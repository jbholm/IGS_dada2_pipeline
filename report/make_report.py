#!/usr/bin/env python3

# -*- coding: utf-8 -*-
import subprocess, glob, argparse, datetime, os, shutil, re, sys, json, urllib.parse, math, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from mako.template import Template
from mako.lookup import TemplateLookup
import htmltag  # Installed manually from github
import PyInquirer  # pip install PyInquirer
import examples
import colour

# ALLOW THIS VARIABLE TO DROP THROUGH INTO CHILD FUNCTIONS
opts = {
    'map': "",
    'tests': ["Per base sequence quality", "Per tile sequence quality",
                  "Sequence Length Distribution"],
        'qcCaptions': {
    "Per base sequence quality": ("", ""),
    "Per tile sequence quality": ("", ""),
    "Sequence Length Distribution": ("Illumina MiSeq machines yield 301-bp reads", "Illumina MiSeq machines yield 301-bp reads")
},
    'contaminants': [],
    'asvDfs': {},  # dataframes containing abundance tables
    'defaultAsvTable': "",
    'controlTable': "",
    'lookup': None
}
    
class Pathwiz(object):
    def __init__(self, scriptsDir, projectDir):
        self.sd = scriptsDir # These are publicly visible
        self.pd = projectDir
    
    # Converts a relative path in the project directory to an abs path
    def proj(self, path):
        return os.path.join(self.pd, path)
    
    # Converts relative path in the script directory to abs path
    def script(self, path):
        return os.path.join(self.sd, path)
    
    # Gets a path relative to the project directory. As a side-effect, makes
    # paths with whitespace safe.
    def relToProj(self, path):
        return str(Path(path).relative_to(self.pd))

class controlsPar(object):
    def __init__(self, title, prefix, patt):
        self.title = title
        self.prefix = prefix
        self.patt = patt
        self.urlsafe = urllib.parse.quote(title)

def initControlsParams():
    ctrlsPars = []
    ctrlsPars.append(controlsPar('PCR Neg. Controls',
                                    "pcrNeg",
                                    "ntcctr|PCRNTC|PCR.NEG|PCR.NTC|PCRNEG|PCRNEGCTRL"))
    ctrlsPars.append(controlsPar('PCR Pos. Controls', "pcrPos",
                                    "PCRPOS|PCR.pos|pCR.pos|POSCTRL|POS.CTRL|POSCON|posctr"))
    ctrlsPars.append(controlsPar('Extraction Neg. Controls', "extNeg",
                                 "EXTNTC|NTC.EXT|NTCEXT|EXT.NTC|NTC|EXTNEG"))
    # Add ext. pos. controls when the time comes
    return(ctrlsPars)
    
def main(pw, args):
    os.chdir(pw.pd)
    reportDir = os.path.join(pw.pd, "REPORT")
    reportFilesDir = os.path.join(reportDir, ".Report_files")
    try:
        overwriteDir(reportDir)
        overwriteDir(reportFilesDir)
    except OSError:
        print(reportDir + " already exists")
    

    opts['lookup'] = TemplateLookup(directories=[scriptDir])
    mytemplate = opts['lookup'].get_template("MSL_REPORT.html")

    map_table = getMapHTML(pw)
    date = datetime.date.today().isoformat()
    dada2stats = getDada2Stats(pw)
    asvTables = getAsvTables(pw)

    # start off the table of contents
    toc = ["""<ul style="background: #ffffff;">
    <li><a href="#info"><span>Sample Information</span></a></li>
            <li><a href="#qc"><span>QC Report</span></a></li>
            <li>Results</li>
            <li><a href="#all-samples"><span class="nav-link-2">All Samples</span></a></li>"""]

    # Get content on controls.
    ctrlPars = initControlsParams() # returns list. Order matters bc the toc is ordered!
    ctrlContent = {}
    for pars in ctrlPars:
        # If this project contained controls, getCtrlPlots returns html and js
        content = getCtrlPlots(pw, pars)
        ctrlContent[pars.title] = content
        if(len(content['js']) > 0):
            toc.append('<li><a href="#' + pars.urlsafe +
                       '"><span class="nav-link-2">' + pars.title + """</span></a></li>""")        

    toc.append("</ul>")
    toc = "\n".join(toc)



    f = open(os.path.join(reportDir, "Report.html"), "w+")
    f.write(mytemplate.render(date=date,
                              mapping=map_table,
                              dada2stats=dada2stats,
                              qc=fastqc(pw),
                              asvtableshtml=asvTables['html'],
                              pcrNegCtrlHtml=ctrlContent['PCR Neg. Controls']['html'],
                              pcrPosCtrlHtml=ctrlContent['PCR Pos. Controls']['html'],
                              extNegCtrlHtml=ctrlContent['Extraction Neg. Controls']['html'],
                              toc=toc))
    f.close()

    mytemplate = opts['lookup'].get_template("data.js")
    
    with open(pw.script("js/data.js"), "w") as f:
        f.write(mytemplate.render(asvtablesjs=asvTables['js'],
                                  pcrNegCtrlJs=ctrlContent['PCR Neg. Controls']['js'],
                                  pcrPosCtrlJs=ctrlContent['PCR Pos. Controls']['js'],
                                  extNegCtrlJs=ctrlContent['Extraction Neg. Controls']['js'],
                contaminants=", ".join([enquote(c) for c in opts['contaminants']])
                )
                )

    myDir = scriptDir

    def cpAssets(srcs):
        for src in srcs:
            srcAbs = os.path.join(myDir, src)
            dest = os.path.join(reportFilesDir, src)
            try:
                shutil.rmtree(dest)
                print("Overwriting " + dest)
            except:
                pass
            shutil.copytree(srcAbs, dest)

    #        files = os.listdir(srcAbs)
    #        for f in files:
    #            shutil.move(os.path.join(srcAbs, f), os.path)
    cpAssets(("js", "css"))

    # minifyCss(pw.pd)

    print("Report created at " + os.path.join(reportDir, "Report.html"))


def overwriteDir(directory):
    try:
        shutil.rmtree(directory)
        print("Overwriting " + directory)
    except OSError:
        pass
    os.makedirs(directory)


def softMkDir(directory):
    try:
        shutil.rmtree(directory)
    except OSError:
        pass

def askFile(dir, msg = "Please choose a file:"):
    oldwd = os.getcwd()
    os.chdir(dir)
    ans = ""
    while len(ans) == 0:
        it = os.scandir()
        choices = list()
        for i, entry in enumerate(os.scandir()):
            letter = "D " if entry.is_dir() else "F "
            choices.append(letter + entry.name)
        choices.sort()
        choices = [".", "..", "Cancel"] + choices
        questions = [
        {
            'type': 'list',
            'name': 'chooser',
            'message': msg,
            'choices': choices,
        }
        ]
        choice = PyInquirer.prompt(questions, style=examples.custom_style_2)['chooser']

        if choice != "." and choice != ".." and choice[0] == "F":
            ans = os.path.join(os.getcwd(), choice[2:len(choice)])
        elif choice == "Cancel":
            ans = ""
            break
        else:
            nav = choice[2:] if choice.startswith("D ") else choice
            os.chdir(nav)
    os.chdir(oldwd)
    return(ans)

def getMappingFile(pw):
    opts['map'] = ""
    oldCwd = os.getcwd()
    os.chdir(pw.pd)
    while len(opts['map']) == 0:
        if args.interactive:
            candidates = glob.glob('**/*mapping*.txt', recursive=True)
            questions = [
            {
                'type': 'list',
                'name': 'useGlobbed',
                'message': "Choose a mapping file",
                'choices': candidates + ["Browse...", "Cancel"]
            }
            ]
            choice = PyInquirer.prompt(questions, style=examples.custom_style_2)['useGlobbed']
            if choice == "Browse...":
                # open file chooser
                opts['map'] = askFile(pw.pd, "Please choose a mapping file")
            elif choice == "Cancel":
                sys.exit(0)
            else:
                opts['map'] = pw.proj(choice)
        else:
            # look recursively in the run directories we were given
            maps = [glob.glob(rundir + '/**/*mapping*.txt', recursive=True)[0] for rundir in args.runs]  
            opts['map'] = pw.proj(joinMaps(maps))
    os.chdir(oldCwd)
    return()

def joinMaps(files): # concatenates files into a tab-delimited mapping file in the current directory
    dest = "project_map.txt"
    with open(dest, "w") as outF:
        with open(files[0], "r") as inF:
            for line in inF:
                outF.write(line.strip() + "\n")
        for map in files[1:len(files)]:
            with open(map, "r") as inF:
                for i, line in enumerate(inF):
                    if i != 0 and re.search(r'^\S+(\t\S+){3}', line) is not None:
                        outF.write(line)
    
    return(dest)

def enquote(s):
    return("\"" + s + "\"")

def getMapHTML(pw):
    status = -1
    while status != 0:
        getMappingFile(pw)
        os.chdir(scriptDir)
        
        cmd = "Rscript " + "mappingToHtml.R -m " + enquote(opts['map'])
        if not args.verbose:
            status = subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if status != 0:
                print("Unable to parse selected file as map.")
                if not args.interactive:
                    sys.exit(1)
                opts['map'] = ""
        else:
            status = subprocess.call(cmd, shell=True)

    f = open("mapping.html.frag", "r")
    map_table = f.read()
    f.close()
    return map_table


def getDada2Stats(pw):
    statsFile = os.path.join(pw.proj("REPORT"), "DADA2_stats.txt")
    oldStatsFile = pw.proj("stats_file_cmp.txt")

    status = -1 # Status of whole operation
    while status != 0:
        status1 = -1
        lines = []
        while status1 != 0:
            try:
                shutil.copy2(oldStatsFile, statsFile)
                with open(oldStatsFile, "r") as f:
                    lines = f.readlines()

                status1 = 0
            except FileNotFoundError:
                print("Couldn't find " + oldStatsFile)
                oldStatsFile = askFile(pw.pd, "Please choose a stats file (combined part1-part2 DADA2 stats)\n")
            except BaseException as e:
                print(str(e))
                oldStatsFile = askFile(pw.pd, "Please choose a stats file (combined part1-part2 DADA2 stats)\n")
        
        status = 0 # Good unless parsing errors encountered in loop
        ans = "\n\t"
        for i in range(len(lines)):
            if i != 0 and status == 0:
                line = lines[i].rstrip()
                fields = line.split("\t")
                if len(fields) == 1:
                    break
                try:
                    ans += "<tr>\n<td>" + fields[0] + "</td>\n<td>" + fields[1] + "</td>\n<td>" + fields[2] + "</td>\n<td>" + \
                        fields[3] + "</td>\n<td>" + fields[4] + "</td>\n</tr>\n"
                except:
                    if args.interactive:
                        print("Parsing DADA2 stats file failed.\nPlease select a different file.")
                        status = 2
                        oldStatsFile = askFile(pw.pd, "Please choose a stats file (combined part1-part2 DADA2 stats)\n")
                        break
                    else:
                        print("Parsing DADA2 stats file failed. Please validate.\n")
                        sys.exit(1)


    ans += "\n"
    return ans


def fastqc(pw):
    ans = ""

    rundirs = ()
    imgDir = os.path.join(pw.proj("REPORT"), "FastQC")
    overwriteDir(imgDir)

    # "rundirs" should have all directories in the project directory
    p = Path(pw.pd)
    rundirs = [pw.proj(run) for run in args.runs]
    # REPORT is a keyword directory that is never considered a run directory
    # Might want to change this so the user can select specific directories as runs.
    runs = 0

    for rundir in rundirs:  # Get FastQC report from each run folder
        p = Path(rundir)

        subdirs = [x for x in p.iterdir() if x.is_dir()]

        # Get the directory names of the forward and reverse FastQC folders.
        fwdFastQc = [os.path.join(x, "seqs_fastqc.zip") for x in subdirs if os.path.basename(x) in ["R1split", "fwdSplit"] and os.path.isfile(os.path.join(x, "seqs_fastqc.zip"))]
        revFastQc = [os.path.join(x, "seqs_fastqc.zip") for x in subdirs if os.path.basename(x) in ["R2split", "R4split", "revSplit"] and os.path.isfile(os.path.join(x, "seqs_fastqc.zip"))]
        # for subdir in subdirs:
        #     
        #     if subdir == "R1split" or subdir == "fwdSplit":
        #         fwdFastQc = os.path.join(subdir, "seqs_fastqc")
        #     elif (subdir in ["R2split", "R4split", "revSplit"]):
        #         revFastQc = os.path.join(subdir, "seqs_fastqc")

        ### UNTESTED ###
        if len(revFastQc) > 1:
            print("Error: multiple FastQC directories for reverse reads in run " + os.path.basename(str(rundir)) + " found.")
            
            questions = [
            {
                'type': 'list',
                'name': 'chooser',
                'message': "Please choose one:",
                'choices': revFastQc,
            }
            ]
            revFastQc = [PyInquirer.prompt(questions, style=examples.custom_style_2)['chooser']]
        if len(fwdFastQc) > 1:
            print("Error: multiple FastQC directories for forward reads in run " + os.path.basename(str(rundir)) + " found.")
            
            questions = [
            {
                'type': 'list',
                'name': 'chooser',
                'message': "Please choose one:",
                'choices': fwdFastQc,
            }
            ]
            fwdFastQc = [PyInquirer.prompt(questions, style=examples.custom_style_2)['chooser']]
        ### Now fwdFastQc and revFastQc are lists containing zero or one filename to a zip folder

        if (len(fwdFastQc) == 1) and (len(revFastQc) == 1):
            runs += 1
            fwdFastQc = str(fwdFastQc[0])
            revFastQc = str(revFastQc[0])

            for zipdir in (fwdFastQc, revFastQc):
                with zipfile.ZipFile(zipdir, 'r') as zip_ref:
                    zip_ref.extractall(os.path.dirname(zipdir[0:len(zipdir) - 4]))
            fwdFastQc = fwdFastQc[0:len(fwdFastQc) - 4]
            revFastQc = revFastQc[0:len(revFastQc)-4]

            def processDir(dir, desc):
                ans = ""
                # COPY all images to the appropriate report folder and get their
                # filepaths (dict of test:filepath pairs)
                imgs = getFastqcLocal(directory=dir)
                # Get the status flags, which have been computed according to our custom rules
                # returns dict of test:status pairs
                statuses = parseFastqcSummary(dir)

                for test in opts['tests']:
                    # Relative to report.html in project directory
                    dirs = imgs[test].split(os.sep)
                    pathRel = os.sep.join(dirs[len(dirs) - 4:len(dirs)])
                    ans += tagFastqc(pathRel) + "\n"
                    ans += captionFastqc(test, statuses[test])
                sectionHead = htmltag.h3(desc)
                paragraph = htmltag.p(htmltag.HTML(ans))
                div = htmltag.div(htmltag.HTML(paragraph))
                return htmltag.div(htmltag.HTML(sectionHead + div), _class="accordion")

            fwd = processDir(fwdFastQc, "Forward reads")
            rev = processDir(revFastQc, "Reverse reads")

            sectionHead = htmltag.h2(os.path.basename(rundir))
            ans += htmltag.div(
                htmltag.HTML(
                    sectionHead +
                    htmltag.div(
                        htmltag.HTML(fwd + rev))),
                _class="accordion accordionRun")
        else:
            if (len(fwdFastQc) == 0):
                checkedDirs = "\n".join([os.path.join(str(rundir), subdir, "seqs_fastqc") for subdir in ["R1split", "fwdSplit"]])
                print("Warning: FastQC files for forward reads in run " + os.path.basename(str(rundir)) + " not found. Looked for:\n" + checkedDirs + "\n")
            if (len(revFastQc) == 0):
                checkedDirs = "\n".join([os.path.join(str(rundir), subdir) for subdir in ["R2split", "R4split", "revSplit"]])
                print("Warning: FastQC files for reverse reads in run " + os.path.basename(str(rundir)) + " not found. Looked for:\n" + checkedDirs + "\n")

    script = htmltag.HTML("""
                          <script>
                    $( ".accordion" ).accordion({collapsible: true, heightStyle: "content"});
                    $( ".accordionRun" ).accordion( "option", "active", 0 );
            </script>""")
    if runs == 1:
        script = script.append("""$( ".accordionRun" ).accordion( "option", "collapsible", false );
                              """)
    return ans + script


'''
COPY all images to the appropriate report folder and get their filepaths
'''


def getFastqcLocal(directory):
    ans = {}
    run = os.path.basename(os.path.dirname(os.path.dirname(directory)))
    imgDir = os.path.join(pw.proj("REPORT"), "FastQC")

    # Get the images in this fastqc analysis
    imgs = [os.path.join(directory, "Images", filename) for filename in
                    ("per_base_quality.png", "per_tile_quality.png",
                     "sequence_length_distribution.png")]

    # R1_fastqc, R2_fastqc, or R4_fastqc, depending on the project directory layout
    dirSpecificDir = os.path.join(imgDir, run, os.path.basename(
        os.path.dirname(directory))[0:2] + "_fastqc")
    os.makedirs(dirSpecificDir)
    if len(imgs) == len(opts['tests']):
        for i in range(len(imgs)):
            file = imgs[i]
            imgDest = os.path.join(dirSpecificDir, os.path.basename(file))
            shutil.copy2(file, imgDest)
            ans[opts['tests'][i]] = imgDest
    return ans


'''
From the fastQC summary.txt, read the test results from four tests:
    Per base sequence quality
    Per tile sequence quality
    Sequence Length Distribution
    Per base N content
'''
def parseFastqcSummary(directory):
    ans = {}

    f = open(os.path.join(directory, "summary.txt"), "r")
    for line in f.readlines():
        fields = line.split("\t")
        for test in opts['tests']:
            if test in fields[1]:
                ans[test] = fields[0]
    f.close
    return ans


'''
For each FastQC image, create an <img> tag that is HTML-ready. Return tags as
a list.
'''
def tagFastqc(img):
    name = os.path.basename(img).split(".")[0]
    altText = " ".join(name.split("_"))
    return htmltag.img(src=img, alt=altText, _class="fastqc")


def captionFastqc(test, status):
    passfail = 1 if status == "PASS" else 0
    caption = opts['qcCaptions'][test][passfail]
    if len(caption) == 0:
        return ""
    if status == "PASS":
        ans = htmltag.div(
            htmltag.HTML(
                htmltag.div(
                    htmltag.HTML(
                        htmltag.p(
                            htmltag.HTML(
                                htmltag.span(
                                    _class="ui-icon ui-icon-info",
                                    style="float: left; margin-right: .3em;"
                                ) +
                                caption
                            )
                        )
                    ),
                    _class="ui-state-highlight ui-corner-all",
                    style="margin-top: 20px; padding: 0 .7em;"
                )),
            _class="ui-widget"
        )
    else:
        ans = htmltag.div(
            htmltag.HTML(
                htmltag.div(
                    htmltag.HTML(
                        htmltag.p(
                            htmltag.HTML(
                                htmltag.span(
                                    _class="ui-icon ui-icon-info",
                                    style="float: left; margin-right: .3em;"
                                ) +
                                caption
                            ))
                    ),
                    _class="ui-state-error ui-corner-all",
                    style="margin-top: 20px; padding: 0 .7em;"
                )),
            _class="ui-widget"
        )
    return ans


def enquote(s):
    return "'" + s + "'"


def readAsvTable(filepath):
    headers = []
    with open(filepath, "r") as f:
        for i in [0, 1]:
            line = f.readline()
            try:
                abundances = line.strip().rstrip(",").split(",")[1:]
                # Next line throws ValueError if header line
                abundances = [int(x) for x in abundances]
            except ValueError:
                headers.append(i)

    # must do this to get columns with duplicate names
    df = pd.read_csv(filepath, header=None, index_col=0) # note header=None
    if df.isnull().sum().sum() > 0:
        print("Warning: Missing values in " + filepath + " have been replaced with 0. Inspect the original CSV with extreme caution!")
        df.fillna(0, inplace = True) # Prevent errors in next step (extra columns should be visible to user in final report)
    
    df.columns = pd.MultiIndex.from_frame(df.iloc[headers].T.astype(str))
    df = df.iloc[ [x not in headers for x in range(len(df.index)) ] ]  # drop header rows
    df = df.astype(int)  # allow dataframe to fix itself
    return(df)


def df_to_js(df, float_format=None):
    items = ["[["]
    for rInd, row in df.iterrows():
        for val in row:
            formatted = float_format.format(
                val) if isinstance(val, float) else val
            items.append(str(formatted) + ",")
        items.append("],\n[")
    items.pop()
    items.append("]]\n")
    return "".join(items)


def getAsvTables(pw):
    os.chdir(pw.pd)
    csvs = glob.glob("*.csv", recursive=True)

    chooseNth = 1
    include = []
    
    if args.interactive:
        while (chooseNth > 0):
    
            if chooseNth == 1:
                question = "Please choose the first table to include"
            else:
                question = "Please choose the next table to include"
        
            questions = [
                {
                    'type': 'list',
                    'name': 'selection',
                    'message': question,
                    'choices': csvs + [PyInquirer.Separator(), "Abort", "Done"]
                }
            ]
        
            selected = PyInquirer.prompt(
                questions, style=examples.custom_style_2)['selection']
        
            if selected == "Abort":
                sys.exit(0)
            elif selected == "Done":
                chooseNth = 0
            else:
                include.append(selected)
                i = csvs.index(selected)
                csvs = csvs[0:i] + csvs[i + 1:]
        
                # Remember which ASV table was chosen first
                if chooseNth == 1:
                    opts['defaultAsvTable'] = selected
        
                chooseNth += 1
                print("\n")
        
    else:
        include = glob.glob(pw.proj('*asvs+taxa.csv'), recursive=True) + glob.glob(pw.proj('*taxa-merged.csv'), recursive=True)
        

     # for selection in include:
    # Get as numeric HTML table if chosen

    # Get as heatmap if chosen

    # f = open(pw.pd + "/IHV_all_runs_dada2_abundance_table_PECAN_taxa_only_merged.csv", "r")
    html = ""
    js = ""
    heatmapNbr = 1
    for file in include:
        fields = os.path.splitext(file)[0].split(sep="_")
        
        active = "active" if heatmapNbr == 1 else ""        
        title = ""
        tbody = ""
        data_str = ""
        taxa = []
        asvIDs = []
        nCol = 0
        xAxes = ""
        data = ""
        zmin = 0
        zmax = 0
        topMargin = 40
        
        shutil.copyfile(os.path.join(file), os.path.join(pw.pd, "REPORT", os.path.basename(file)))

        taxnmy = fields[len(fields)-2]
        if taxnmy == "PECAN-SILVA":
            taxnmy += "*"
        
        form = fields[len(fields)-1]
        if form in ("taxa", "taxa+asvs", "asvs+taxa"):
            title="ASV counts labeled with " + taxnmy + " taxonomic assignments"
        elif form == "taxa-merged":
            title = "Counts by " + taxnmy + " taxonomic assignment"
        else:
            print("Error: can't determine count table file format from filename.")
        
        opts['asvDfs'][file] = readAsvTable(pw.proj(file))


        # Get the one or two headers in the dataframe and determine whether
        # they are ASV IDs or taxa
        for i in range(2):
            try:
                colnames = list(
                    opts['asvDfs'][file].columns.get_level_values(i))
                if all([str(name)[0:3] == "ASV" for name in colnames]):
                    asvIDs = colnames
                else:
                    taxa = colnames
            except IndexError:
                break
        sampleIDs = [enquote(id) for id in list(opts['asvDfs'][file].index)]

        # Create an HTML table body for displaying a table of samples, total reads, and top hits
        tdata = []
        for index, row in opts['asvDfs'][file].iterrows():
            tophiti = row.values.argmax() # index of top hit

            if isinstance(tophiti, slice):
                tophiti = str(tophiti)
                tophiti = int(re.search("(?<=slice\()[0-9]+", str(tophiti)).group(0))
            topASV = asvIDs[tophiti] if len(asvIDs) > 0 else ""
            topTaxon = taxa[tophiti] if len(taxa) > 0 else ""

            tdata.append([index, sum(row), [topASV, topTaxon]])

        def maketd(data):
            return ("""<td style="text-align:left;">""" + str(data) + "</td>\n")

        def makerow(row, i):
            tr = ["<tr>"]
            tr.append(maketd(row[0]))
            tr.append(maketd(row[1]))
            sample = row[0]
            asv = row[2][0]
            taxon = row[2][1]
            taxonStr = " (" + taxon + ")" if(len(taxon) > 0) else ""
            tr.append(maketd("<a onclick=\"asvHeatmapMaxHighlight(" + str(i) + ", 'heatmapInner-" + str(heatmapNbr) + "')\" href=\"#\">" + asv + taxonStr + "</a>"))
            tr.append("</tr>")
            return("\n".join(tr))
        tbody = "\n".join([makerow(row, i) for i, row in enumerate(tdata)])

        # Add quotes to each xlabel (for javascript syntax)
        taxa = [enquote(taxon) for taxon in taxa]
        asvIDs = [enquote(asvID) for asvID in asvIDs]

        #data_np = np.array(data_arr, dtype=float)
        rowsums = opts['asvDfs'][file].sum(axis=1)
        relAbund = opts['asvDfs'][file].div(rowsums, axis="index")
        zmin = relAbund.min().min()
        zmax = relAbund.max().max()
        zmax = zmax if zmax < 0.3 else np.ceil(zmax * 10) / 10
        
        # Reorder by sum(relative abundance) across samples
        colsums = relAbund.sum(axis=0)
        relAbundOrder = list(colsums.argsort().iloc[::-1])
        relAbund = relAbund.iloc[:, relAbundOrder]
        if len(taxa) > 0:
            taxa = [taxa[i] for i in relAbundOrder]
            
            
        if len(asvIDs) > 0:
            asvIDs = [asvIDs[i] for i in relAbundOrder] 
            topMargin = 100
        
        data_str = df_to_js(
            relAbund, float_format='{:.' + str(int(np.amax(np.ceil(np.log10(rowsums))))) + 'f}')
        #data_str = "[[" + data_str[0:(len(data_str) - 4)] + "]]\n"
        nCol = len(taxa) if len(taxa) > 0 else len(asvIDs)
        
        
        
        relAbundOrder = "[" + ", ".join([str(x) for x in relAbundOrder]) + "]\n"

        htmlTemplate = opts['lookup'].get_template("countTable.html")
        jsTemplate = opts['lookup'].get_template("countTable.js")
        data = "["
        if len(asvIDs) > 0:
            xAxes += """xaxis2: {
                ticktext: xTopValues,
                tickvals: xTickPos,
                tickfont: {color: "#333"},
                title: '<b>ASVs</b>',
                side: 'top',
                type: 'category',
                automargin: true,
                fixedrange: true,
                linecolor: '#333',
    linewidth: 1,
            },"""
            data += """{  // second trace for top x-axis
                       x: xTopValues,
                       y: yValues,
                       z: zValues,
                       type: 'heatmap',
                       mode: 'markers',
                       colorscale: 'Viridis',
                        showscale: false,
                       xaxis: 'x2',
                       hoverinfo: "text",
                   text: zValues.map((row, i) => row.map((item, j) => {
                       return "Sample: " + yValues[i] + "<br>" + xTopValues[j] + "<br>Taxon: " + xBotValues[j] + "<br>Rel. abundance: " + item;
                       }))
                   },"""
        if len(taxa) > 0:
            xAxes += """xaxis: {
                ticktext: xBotValues,
                tickvals: xTickPos,
                tickfont: {color: "#333"},
                title: "<b>Assigned taxa</b>",
                showgrid: false,
                type: 'category',
                automargin: true,
                fixedrange: true,
                linecolor: '#333',
    linewidth: 1,
            },"""
            data += """{
            x: xBotValues,
            y: yValues,
            z: zValues,
            type: 'heatmap',
            mode: 'markers',
            colorscale: 'Viridis',
            showscale: false,
            hoverinfo: "text",
        text: zValues.map((row, i) => row.map((item, j) => {
            return "Sample: " + yValues[i] + "<br>Taxon: " + xBotValues[j] + "<br>Rel. abundance: " + item;
            }))
        },"""
        if len(taxa) == 0 and len(asvIDs) == 0:
            print("WARNING: Abundance tables did not have headers. Plotting might be broken.")
        data += "]"

        html += htmlTemplate.render(active=active,
                                    hmNbr=heatmapNbr,
                                    title=title,
                                    tbody=tbody,
                                    isAsvTable=len(asvIDs) > 0) + "\n"
        js += jsTemplate.render(abundances=data_str,
                                taxa=", ".join(taxa), asvs=", ".join(asvIDs),
                                samples=", ".join(sampleIDs),
                                nCol=nCol,
                                xAxes=xAxes, data=data,
                                zmin=zmin, zmax=zmax,
                                height = (15 * len(sampleIDs) + 250 + topMargin),
                                topmargin = topMargin,
                                keyTopMargin = topMargin - 20,
                                hmNbr=heatmapNbr)

        heatmapNbr += 1
        
    sectionTemplate = opts['lookup'].get_template("countPage.html")
    selectorOpts = ""
    for id, filepath in enumerate(include):
        selectorOpts += '<option value="heatmapTab-' + \
            str(id + 1) + '">...' + pw.relToProj(filepath) + '</option>\n'
    
    sectionHtml = sectionTemplate.render(tables=html, selectOptions=selectorOpts)

    return {'html': sectionHtml, 'js': js}

def getCtrlPlots(pw, pars):
    if len(opts['asvDfs']) == 0:
        return {'html': "", 'js': ""}

    if(len(opts['controlTable']) == 0):
        if args.interactive:
            question = "Extract control data from " + \
                opts['defaultAsvTable'] + " ?"
    
            questions = [
                {
                    'type': 'confirm',
                    'name': 'useDefault',
                    'message': question,
                    'default': True,
                }
            ]
    
            if PyInquirer.prompt(questions, style=examples.custom_style_2)['useDefault']:
                df = opts['defaultAsvTable']
            else:
                question = "Please choose the ASV table from which to extract control data:"
                questions = [
                    {
                        'type': 'list',
                        'name': 'selection',
                        'message': question,
                        'choices': opts['asvDfs'].keys()
                    }
                ]
                df = PyInquirer.prompt(questions, style=examples.custom_style_2)[
                    'selection']
    
            opts['controlTable'] = df
        else:
            opts['controlTable'] = glob.glob(pw.proj('*asvs+taxa.csv'), recursive=True)[0]


    # Get control samples
    mapping = pd.read_csv(opts['map'],
                          sep="\t", header=0, index_col=3)
    ctrls = list(mapping.filter(regex=pars.patt, axis=0).iloc[:, 0])
    if(len(ctrls) == 0):
        return {'html': '',
                'js': ''}

    # Get the taxa/ASV IDs
    asvIDs = []
    taxa = []
    for i in range(2):
        try:
            colnames = list(
                opts['asvDfs'][opts['controlTable']].columns.get_level_values(i))
            if all([name[0:3] == "ASV" for name in colnames]):
                asvIDs = colnames
            else:
                taxa = colnames
        except IndexError:
            break
    if(len(asvIDs) == 0):
        print("Warning: " + opts['controlTable'] + \
            " does not have ASV IDs. Unable to highlight ASVs appearing in " + \
            template_prefix + " controls.")

    # Get their data
    data = opts['asvDfs'][opts['controlTable']].filter(items=ctrls, axis=0)
    meds = list(data.median(axis=0))
    sems = list(data.sem(axis=0).div(math.sqrt(len(ctrls))))

    # subset the taxa with mean abundance > 0
    idxs = [i for i, x in enumerate(meds) if x > 0]
    meds = [meds[i] for i in idxs]
    data = data.iloc[:, idxs]
    if(len(taxa) > 0):
        taxa = [taxa[i] for i in idxs]
    sems = [sems[i] for i in idxs]
    if(len(asvIDs) > 0):
        asvIDs = [asvIDs[i] for i in idxs]

    # (now sort by mean abundance)
    abundOrder = list(np.argsort(meds))
    abundOrder.reverse()
    if(len(taxa) > 0):
        taxa = [taxa[i] for i in abundOrder]
    sems = [sems[i] for i in abundOrder]
    data = data.iloc[:, abundOrder]
    if(len(asvIDs) > 0):
        asvIDs = [asvIDs[i] for i in abundOrder]

    opts['contaminants'] = asvIDs

    html_template =opts['lookup'].get_template(pars.prefix + "CtrlPlot.html")
    sample_names = data.index.values.tolist()
    options = ['<option>' + name + '</option>' for name in sample_names]
    colors = list(colour.Color("#F0E442").range_to(
        colour.Color("#CC79A7"), len(sample_names)))
    colors = ",".join([enquote(col.hex) for col in colors])

    template = opts['lookup'].get_template(pars.prefix + "CtrlPlot.js")
    js = template.render(xs=", ".join([str(x) for x in range(len(taxa))]),
                         ys=df_to_js(data.T),
                         taxa=", ".join([enquote(taxon) for taxon in taxa]),
                         asvs=", ".join([enquote(asvID) for asvID in asvIDs]),
                         sems=", ".join(
                             ["{:.2f}".format(x) for x in sems]),
                         ntaxa=str(len(taxa)), colors=colors)

    html = html_template.render(id=pars.urlsafe,
                                options="\n".join(options))
    return {'html': html, 'js': js}


def minifyCss(pw):
    #pw = pw.proj("REPORT")
    #subprocess.check_call("css-purge -i " + enquote(pw.proj(".Report_files", "css")) + " -m " +
                          #enquote(pw.proj("Report.html")) + " -o " + enquote(pw.proj("report.css")), shell=True)
    #shutil.rmtree(pw.proj(".Report_files", "css"))
    #os.mkdir(pw.proj(".Report_files", "css"))
    #shutil.move(pw.proj("report.css"), os.path.join(
        #pw, ".Report_files", "css", "report.css"))

    #shutil.move(pw.proj("Report.html"),
                #pw.proj("Report.html.bak"))
    #with open(pw.proj("Report.html.bak")) as inFile:
        #with open(pw.proj("Report.html"), "w+") as outFile:
            #patt = re.compile(
                #'<link[^<]*(?=rel="stylesheet")[^>]*>', flags=re.DOTALL)
            #contents = inFile.read()
            #contents = re.sub(patt, "", contents, count=0)
            #patt = re.compile('</head>')
            #outFile.write(re.sub(patt, string=contents,
                                 #repl="\n<link href='.Report_files/css/report.css' rel='stylesheet'>\n</head>"))

    #os.remove(pw.proj("Report.html.bak"))
    return


if __name__ == '__main__':
    scriptDir = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(
        description='''A script that automatically populates fields in the HTML report template. ''',
        epilog="""""")
    parser.add_argument("wd", nargs=1, metavar="PROJECT_DIR",
                        help="The project directory containing a map, stats_file_cmp.txt, any ASV tables in CSV format, and run folder(s) with FastQC reports.")
    parser.add_argument("--verbose", action='store_true')
    parser.add_argument("--interactive", "-i", action='store_true')
    parser.add_argument("--runs", "-r", nargs='+', help="The names of the run directories to include. These must be subdirectories of DIR. If not --interactive, --run is required.")
    args = parser.parse_args()
    np.set_printoptions(threshold=np.inf)
    if not args.interactive and args.runs is None:
        print("Requires either --interactive or --runs")
        sys.exit(1)
    for wd in args.wd:
        pw = Pathwiz(scriptsDir=scriptDir,
                     projectDir=wd)
        main(pw, args=args)
