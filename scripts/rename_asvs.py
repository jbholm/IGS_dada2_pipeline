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
parser.add_argument('--pecan', nargs="?", const="", default="",
                    help='A pecan classification file (usu. *MC_order7_results.txt)')
parser.add_argument('--classif', '-c', nargs='*', help='Any non-PECAN classification CSVs containing one header line, ASVs in the first column and taxonomic assignments in the remaining columns. Do NOT omit the project prefix from these filenames, if present.')
args = parser.parse_args()

csv = args.project + "_all_runs_dada2_abundance_table.csv"
fasta = "all_runs_dada2_ASV.fasta"
taxaCsv = args.project + "_all_runs_dada2_abundance_table_w_taxa.csv"
if args.classif is None:
    args.classif = []
classifs = []
for file in args.classif:
    classifs.append(file)

# There are three files containing headers/rownames that reference sequences
# in a FASTA. Are the sequences, headers, and rownames in the same order?

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

csTaxaCSV = ""
with open(taxaCsv, "r") as file:
    # comma-separated string of asvs from the CSV
    csTaxaCSV += file.readline().strip()

csTaxaCSV = csTaxaCSV[1:(len(csTaxaCSV) - 1)]
# Remove trailing and leading commas

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


pecanIdent = True
if args.pecan != "":
    with open(args.pecan, "r") as file:
        csPecan = ""
        for line in file:
            csPecan += line.split()[0] + ","

        csPecan = csPecan[0:(len(csPecan) - 1)]
    if csPecan != csCSV:
        pecanIdent = False

# Test if all the classification files have ordered ASVs identical to the
# abundance table
classifsIdent = True
for csClassif in csClassifs:
    if csClassif != csCSV:
        classifsIdent = False
# Test that the abundance tables and reference FASTA have identical ordered ASVs
if any([csCSV != csFASTA, csTaxaCSV != csCSV, not classifsIdent, not pecanIdent]):
    print >> sys.stderr, "ASVs in the DADA2 unclassified ASV abundance table and the ASV fasta are not in the same order (or might not even be identical sets of ASVs).\nModify scripts/rename_asvs.py to handle this scenario.\n\n"
    sys.exit(1)

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

shutil.move(taxaCsv, taxaCsv + ".bak")
baks.append(taxaCsv + ".bak")
try:
    with open(taxaCsv + ".bak", "r") as inFh:
        inFh.readline()  # discard old header
        # for some reason, this file has an extra column on the header
        csAsvIds = "," + ",".join(asvIds) + ",\n"

        # Write to new csv, and then copy from the backup to the new file
        with open(taxaCsv, "w") as outFh:  # write-after-truncate mode
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
                            ",".join(fields[1:len(fields) - 1]) + "\n"
                    outFh.write(line)
                    ln += 1
    except Exception as error:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        # Restore everything to the way it was
        revertBaks(baks)
        print str(exc_tb.tb_lineno) + ": " + str(error)
        sys.exit(1)

if args.pecan != "":
    shutil.move(args.pecan, args.pecan + ".bak")
    baks.append(args.pecan + ".bak")
    try:
        with open(args.pecan + ".bak", "r") as inFH:
            with open(args.pecan, "w") as outFH:
                ln = 0
                for line in inFH:
                    fields = line.split()
                    fields[0] = asvIds[ln]
                    outFH.write("\t".join(fields) + "\n")
                    ln += 1
    except Exception as error:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        # Restore everything to the way it was
        revertBaks(baks)
        print str(exc_tb.tb_lineno) + ": " + str(error)
        sys.exit(1)

removeBaks(baks)
