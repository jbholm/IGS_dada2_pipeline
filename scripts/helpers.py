import subprocess, argparse, os, sys, re

def readTwoColumnASVKey(file, header = None, phylocolumn = 1){
    asvTaxa = {}
    lineno = 0
    with open(file, "r") as inFile:
        line = inFile.readline()
        if lineno == 0:
            # for the first line, if header wasn't specified (or specified as None), figure out if header
            # this requires reading in two lines
            secondLine = inFile.readline()
            fields1 = line.strip().split()
            fields2 = secondLine.strip().split()
            if header is None:
            
                if len(fields1) == len(fields2):
                    header = False
                else:
                    header = True
            if not header:
                fields = fields1[0:2]
                asvTaxa{fields[0]} = fields[phylocolumn]
                # if no header, do what's necessary to parse the first line anyway
            line = secondLine            

        while (len(line) > 0):
            line = line.strip()
            fields = line.split()[0:2]
            asvTaxa{fields[0]} = fields[phylocolumn]
            line = inFile.readline()
    
    return(asvTaxa)
}

def readFullASVKey(file){
    asvTaxa = {}

    with open(file, "r") as inFile:
        line = inFile.readline()
        
        levelHints = line.strip().split(",")
        speciesCols = [col if re.search(r'(?i)^species$', ) for col in range(len(levelHints))]
        genusCols = [col if re.search(r'(?i)^genus$', ) for col in range(len(levelHints))]
        levelHints = [level[0] + "_" for level in levelHints]

        line = inFile.readline()
        while (len(line) > 0):
            line = line.strip()
            fields = line.split()

            if len(speciesCols) > 0: 
                taxon = fields[speciesCols[0]]
                if len(genusCols) > 0:
                    genus = fields[genusCols[0]]
                    taxon = genus + "_" + species
            elif len(fields) > 1:
                taxon = fields[len(fields) - 1]
                taxon = levelHints[len(fields) - 1] + "_" + taxon
            else:
                taxon = "unclassified"
            asvTaxa{fields[0]} = taxon
        
    return(asvTaxa)
}
