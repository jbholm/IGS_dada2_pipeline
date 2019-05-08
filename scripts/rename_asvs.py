#!/usr/bin/env python2
import shutil, sys, math

csv = sys.argv[1] + "_all_runs_dada2_abundance_table.csv"
fasta = "all_runs_dada2_ASV.fasta"
taxaCsv = sys.argv[1] + "_all_runs_dada2_abundance_table_w_taxa.csv"
classif = sys.argv[1] + "_silva_classification.csv"

# There are three files containing headers/rownames that reference sequences 
# in a FASTA. Are the sequences, headers, and rownames in the same order?

file = open(csv, "r")
csCSV = file.readline().strip() # comma-separated string of asvs from the CSV
file.close()

csFASTA = ""
file = open(fasta, "r")
ln = 1
for line in file:
    if ( ln % 2 == 0 ):
        line = line.strip()
        csFASTA += "," + line # Join all even lines by a comma
    ln += 1
file.close()

csTaxaCSV = ""
with open(taxaCsv, "r") as file:
    csTaxaCSV += file.readline().strip() # comma-separated string of asvs from the CSV
    csTaxaCSV = csTaxaCSV[0:(len(csTaxaCSV) - 1)]

csClassif = ""
with open(classif, "r") as file:
    ln = 1
    for line in file:
        if ln != 1:
            csClassif += "," + line.strip().split(",")[0]
        ln += 1

if any([csCSV != csFASTA, csTaxaCSV != csCSV, csClassif != csCSV]):
    print >> sys.stderr, "ASVs in the DADA2 unclassified ASV abundance table and the ASV fasta are not in the same order (or might not even be identical sets of ASVs).\nModify scripts/rename_asvs.py to handle this scenario.\n\n"
    sys.exit(1)

# Create a backup of the file so we can change just the first line
shutil.move(csv, csv + ".bak")
queue = ["shutil.move(csv + \".bak\", csv)"]
def execute_queue(queue):
    for cmd in queue:
        eval(cmd)

try:
    asvIds = []
    with open(csv + ".bak", "r") as inFh:
        csCSV = inFh.readline().strip()  # OLD comma-separated string of asvs from the CSV

        # Count the ASVs so they can be given _sortable_ IDs
        asvCount = len(csCSV.split(",")) - 1
        for i in range(asvCount):
            asvIds.append("ASV" + str(i).zfill(int(math.ceil(math.log10(asvCount)))))
        csAsvIds = "," + ",".join(asvIds) + "\n"

        # Write to new csv, and then copy from the backup to the new file
        with open(csv, "w") as outFh: # write-after-truncate mode
            outFh.write(csAsvIds)
            line = inFh.readline()
            while (len(line) > 0):
                outFh.write(line)
                line = inFh.readline()
except Exception as error:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    execute_queue(queue)
    print str(exc_tb.tb_lineno) + ": " + str(error)
    sys.exit(1)

shutil.move(taxaCsv, taxaCsv + ".bak")
queue.append("shutil.move(taxaCsv + \".bak\", taxaCsv)")
try:
    with open(taxaCsv + ".bak", "r") as inFh:
        inFh.readline() # discard old header
        csAsvIds = "," + ",".join(asvIds) + ",\n" # for some reason, this file has an extra column on the header

        # Write to new csv, and then copy from the backup to the new file
        with open(taxaCsv, "w") as outFh: # write-after-truncate mode
            outFh.write(csAsvIds)
            line = inFh.readline()
            while (len(line) > 0):
                outFh.write(line)
                line = inFh.readline()
except Exception as error:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    execute_queue(queue)
    print str(exc_tb.tb_lineno) + ": " + str(error)
    sys.exit(1)

shutil.move(fasta, fasta + ".bak")
queue.append("shutil.move(fasta + \".bak\", fasta)")
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
    execute_queue(queue)
    print str(exc_tb.tb_lineno) + ": " + str(error)
    sys.exit(1)

shutil.move(classif, classif + ".bak")
queue.append("shutil.move(classif + \".bak\", classif)")
try:
    with open(classif + ".bak", "r") as inFh:
        with open(classif, "w") as outFh:
            ln = 1
            for line in inFh:
                 # Copy header verbatim
                 # For every other line, substitute in the new ASV ID
                if ln != 1:
                    fields = line.split(",")
                    line = asvIds[ln - 2] + "," + ",".join(fields[1:len(fields) - 1]) + "\n"
                outFh.write(line)
                ln += 1
except Exception as error:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    # Restore everything to the way it was
    execute_queue(queue)
    print str(exc_tb.tb_lineno) + ": " + str(error)
    sys.exit(1)



