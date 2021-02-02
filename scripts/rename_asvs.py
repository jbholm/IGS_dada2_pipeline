#!/usr/bin/env python2
import shutil
import sys
import math
import os
import argparse

parser = argparse.ArgumentParser(
    description='''A script that renames ASVs in a 16S project. ASVs are given unique IDs in the template ASV###, where ### is a zero-padded number of sufficient length to uniquely identify all ASVs.''',
    epilog="""""")
parser.add_argument('--project', '-p', required=True,
                    help='The project name (the prefix for the abundance tables)')
parser.add_argument('--twocol', nargs="*", metavar = "*.txt", help='Any two-column classification files (usu. *MC_order7_results.txt)')
parser.add_argument('--classif', '-c', nargs='*', metavar = "*.csv", help='Any non-PECAN classification CSVs containing one header line, ASVs in the first column and taxonomic assignments in the remaining columns. Do NOT omit the project prefix from these filenames, if present.')
parser.add_argument('fasta', metavar="*.fasta", help='A FASTA file')
parser.add_argument('--no-match-seqs', action='store_true')
args = parser.parse_args()

csv = args.project + "_all_runs_dada2_abundance_table.csv"
if args.fasta:
    fasta = args.fasta
#taxaCsv = args.labeledCsv
classifs = []
if args.classif:
    for file in args.classif:
        classifs.append(file)
twoCols = []
if args.twocol:
    for file in args.twocol:
        twoCols.append(file)

# There are three files containing headers/rownames that reference sequences
# in a FASTA. Are the sequences, headers, and rownames in the same order?

if not args.no_match_seqs:

    file = open(csv, "r")
    csCSV = file.readline().strip()  # comma-separated string of asvs from the CSV
    csCSV = csCSV[1:len(csCSV)]  # Remove a comma at the beginning of the line
    file.close()

    csFASTA = ""
    file = open(fasta, "r")
    ln = 1
    for line in file:
        if (ln % 2 == 1):
            line = line.strip()
            line = line[1:len(line)]
            csFASTA += "," + line  # Join all even lines by a comma
        ln += 1
    file.close()
    # Remove a comma at the beginning of the line
    csFASTA = csFASTA[1:len(csFASTA)]

    # csTaxaCSV = ""
    # with open(taxaCsv, "r") as file:
    #     csTaxaCSV += file.readline().strip() # comma-separated string of asvs from the CSV

    # csTaxaCSV = csTaxaCSV[1:(len(csTaxaCSV))]
    # Remove leading comma

    csClassifs = []
    for classif in classifs:
        csClassif = ""
        with open(classif, "r") as file:
            ln = 1
            for line in file:
                if ln != 1:
                    csClassif += "," + line.strip().split(",")[0]
                ln += 1
        csClassif = csClassif[1:len(csClassif)]
        csClassifs.append(csClassif)

    csTwoCols = []
    for twoCol in twoCols:
        csTwoCol = ""
        with open(twoCol, "r") as file:
            for line in file:
                csTwoCol += line.split()[0] + ","

            csTwoCol = csTwoCol[0:(len(csTwoCol) - 1)]
        csTwoCols.append(csTwoCol)
        
        # for asvPecan, asvCSV in zip(csPecan.split(","), csCSV.split(",")):
        #     if asvPecan != asvCSV:
        #         pecanIdent = False
        #         print(asvPecan + " in " + args.pecan + " is in the same position as " + asvCSV + " in " + csv + "!")


    # Test if all the classification files have ordered ASVs identical to the
    # abundance table
    twoColsIdent = True
    for csTwoCol in csTwoCols:
        for asvClassif, asvCSV in zip(csTwoCol.split(","), csFASTA.split(",")):
            if (asvClassif.strip() != asvCSV.strip()):
                msg = asvClassif + " in " + ",".join(twoCols) + " is in the same position as " + asvCSV + " in " + csv + "!"
                print >> sys.stderr, msg
                twoColsIdent = False

    classifsIdent = True
    for csClassif in csClassifs:
        for asvClassif, asvCSV in zip(csClassif.split(","), csFASTA.split(",")):
            if (asvClassif.strip() != asvCSV.strip()):
                msg = asvClassif + " in " + ",".join(classifs) + " is in the same position as " + asvCSV + " in " + csv + "!"
                print >> sys.stderr, msg
                classifsIdent = False
    csvIdent = True
    for asvCSV, asvFasta in zip(csCSV.split(","), csFASTA.split(",")):
        if (asvCSV.strip() != asvFasta.strip()):
            msg = asvCSV + " in " + csv + " is in the same position as " + asvFasta + " in " + fasta + "!"
            print >> sys.stderr, msg
            csvIdent = False
    # Test that all files provided with command-line options match the provided fasta
    if any([not csvIdent, not classifsIdent, not twoColsIdent]):
        print >> sys.stderr, "ASVs in the DADA2 unclassified ASV abundance table and the ASV fasta are not in the same order (or might not even be identical sets of ASVs).\nModify scripts/rename_asvs.py to handle this scenario.\n\n"
        # sys.exit(1)

# Create a backup of the file so we can change just the first line
shutil.move(csv, csv + ".bak")
# Start a baks that we can use to rollback all changes or delete all backups
baks = [csv + ".bak"]


def revertBaks(baks):
    for file in baks:
        shutil.move(file, file[0:len(file) - 4])


def removeBaks(baks):
    for file in baks:
        os.remove(file)


asvIds = []
try:
    with open(csv + ".bak", "r") as inFh:
        csCSV = inFh.readline().strip()  # OLD comma-separated string of asvs from the CSV

        # Count the ASVs so they can be given _sortable_ IDs
        asvCount = len(csCSV.split(",")) - 1
        for i in range(asvCount):
            asvIds.append(
                "ASV" + str(i).zfill(int(math.ceil(math.log10(asvCount)))))
        csAsvIds = "," + ",".join(asvIds) + "\n"

        # Write to new csv, and then copy from the backup to the new file
        with open(csv, "w") as outFh:  # write-after-truncate mode
            outFh.write(csAsvIds)
            line = inFh.readline()
            while (len(line) > 0):
                outFh.write(line)
                line = inFh.readline()
except Exception as error:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    revertBaks(baks)
    print str(exc_tb.tb_lineno) + ": " + str(error)
    sys.exit(1)

# shutil.move(taxaCsv, taxaCsv + ".bak")
# baks.append(taxaCsv + ".bak")
# try:
#     with open(taxaCsv + ".bak", "r") as inFh:
#         inFh.readline() # discard old header
#         csAsvIds = "," + ",".join(asvIds) + ",\n" # for some reason, this file has an extra column on the header

#         # Write to new csv, and then copy from the backup to the new file
#         with open(taxaCsv, "w") as outFh: # write-after-truncate mode
#             outFh.write(csAsvIds)
#             line = inFh.readline()
#             while (len(line) > 0):
#                 outFh.write(line)
#                 line = inFh.readline()
# except Exception as error:
#     exc_type, exc_obj, exc_tb = sys.exc_info()
#     revertBaks(baks)
#     print str(exc_tb.tb_lineno) + ": " + str(error)
#     sys.exit(1)

shutil.move(fasta, fasta + ".bak")
baks.append(fasta + ".bak")
try:
    with open(fasta + ".bak", "r") as inFh:
        with open(fasta, "w") as outFh:
            ln = 1
            for line in inFh:
                if (ln % 2 == 1):  # odd-number lines contain fasta headers
                    # We know the FASTA has ASVs in the same order as the CSV
                    line = ">" + asvIds[ln / 2] + "\n"
                outFh.write(line)
                ln += 1
except Exception as error:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    # Restore everything to the way it was
    revertBaks(baks)
    print str(exc_tb.tb_lineno) + ": " + str(error)
    sys.exit(1)

for classif in classifs:
    shutil.move(classif, classif + ".bak")
    baks.append(classif + ".bak")
    try:
        with open(classif + ".bak", "r") as inFh:
            with open(classif, "w") as outFh:
                ln = 1
                for line in inFh:
                    # Copy header verbatim
                    # For every other line, substitute in the new ASV ID
                    if ln != 1:
                        fields = line.split(",")
                        line = asvIds[ln - 2] + "," + \
                            ",".join(fields[1:len(fields)])
                    outFh.write(line)
                    ln += 1
    except Exception as error:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        # Restore everything to the way it was
        revertBaks(baks)
        print str(exc_tb.tb_lineno) + ": " + str(error)
        sys.exit(1)

for twoCol in twoCols:
    shutil.move(twoCol, twoCol + ".bak")
    baks.append(twoCol + ".bak")
    try:
        with open(twoCol + ".bak", "r") as inFh:
            with open(twoCol, "w") as outFh:
                ln = 0
                for line in inFh:
                
                    fields = line.split()
                    fields[0] = asvIds[ln]
                    outFh.write("\t".join(fields) + "\n")
                    ln += 1
    except Exception as error:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        # Restore everything to the way it was
        revertBaks(baks)
        print str(exc_tb.tb_lineno) + ": " + str(error)
        sys.exit(1)

removeBaks(baks)
