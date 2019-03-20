# IGS_dada2_pipeline
Pipeline for processing 16S rRNA or ITS amplicon sequences specifically for the University of Maryland Baltimore Institute for Genome Sciences

illumina_dada2.pl is Part 1 for processing samples through the IGS 16S pipeline. 

THE FOLLOWING STEPS ARE REQUIRED BEFORE THIS SCRIPT WILL RUN:

  1. Start from a Rhel7 machine OR Qlogin to RHEL7:
        qlogin -P jravel-lab -l mem_free=500M -q interactive.q
  2. Then enter:
        export LD_LIBRARY_PATH=/usr/local/packages/gcc/lib64
        source /usr/local/packages/usepackage/share/usepackage/use.bsh
        use python-2.7
        use qiime
  3. then run script as below:
  
FOR 2-STEP:
If raw reads are labeled as R1, R2, R3, and R4:
  illumina_dada2.pl -i <directory containing raw reads> -p <project name> -r <run ID> -m <mapping file> -v <variable region> -sd <storage directory either "scratch" or "groupshare" (don't include quotes)>
  
  OR
  
If raw reads are labeled as R1, i1, i2, and R2:
  illumina_dada2.pl -r1 <full path to raw R1 file> -r2 <full path to raw i1 file>
                    -r3 <full path to raw i2 file> -r4 <full path to raw R2 file>
                    -p <project name> -r <run ID> -m <mapping file> -v <variable region>
                    -sd <storage directory>
  
FOR 1-STEP
  add --1Step

For qsub:
  qsub -cwd -b y -l mem_free=1G -P jravel-lab -q threaded.q -pe thread 4 -V -e <path_to_logs>
  -o <path_to_logs>  illumina_dada2.pl -i <directory containing raw reads> -p <project name> -r <run ID> -m <mapping file> -v <variable region> -sd <storage directory enter the word scratch or groupshare>
  
part2_illumina_dada2.pl is Part 2 for processing samples through the IGS 16S Pipeline
  1. Start from a Rhel7 machine OR Qlogin to RHEL7:
        qlogin -P jravel-lab -l mem_free=500M -q interactive.q
  2. Then enter:
        export LD_LIBRARY_PATH=/usr/local/packages/gcc/lib64
        source /usr/local/packages/usepackage/share/usepackage/use.bsh
        use python-2.7
        use qiime
  3. then run script as below from within the project directory created in Step 1:
part2_illumina_dada2.pl -i <comma-separated-input-run-names> -v <variable-region> -p <project-ID>

# To use an earlier version of this pipeline:
## Desktop
Find the desired version on https://github.com/jbholm/IGS_dada2_pipeline/commits/master. Click the "**<>**" button to browse the repository. Click the green "**Clone or download**" dropdown and click "**Download ZIP**" to download this version of the pipeline.
## Terminal
Navigate to a local clone of this repository (should be a directory named "IGS_dada2_pipeline" containing a subdirectory ./.Git/). Enter the following commands to remove all local un-commited changes and checkout a previous commit (change **1.0.0** to your desired version).  
`git add .`  
`git reset --hard`  
`git rev-list feature/version --grep='v1.0.0' | head -n 1 | xargs git checkout`  
If you do not have a local clone already, you will need to provide your GitHub credentials to clone the repository at a command prompt. Follow the desktop instructions to browse the desired version of the repository. Click the green "**Clone or download**" dropdown and copy the HTTPS URL. Then run:
`git clone PASTE/HTTPS/URL/HERE`