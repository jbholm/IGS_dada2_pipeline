# IGS_dada2_pipeline

Two-part pipeline for processing 16S rRNA or ITS amplicon sequences specifically for the University of Maryland Baltimore Institute for Genome Sciences


## Pre-processing 16S or ITS paired-end reads:
`part1.sh` demultiplexes, filters, trims, and merges 16S or ITS paired-end reads. The script can be launched from any location on the IGS server according to this template:

`part1.sh (-i <input directory> | -r1 <fwd reads> -r2 <rev reads> [-i1 <index 1> -i2 <index 2>]) -p <project> -r <run> -m <map> [-v <variable region>] [--1Step] [<options>]`

Please run `part1.sh --help` for complete documentation. 

Running `part1.sh` outside of IGS servers is not yet supported. It is possible however to modify hard-coded paths in `part1.sh` and `illumina_dada2.pl` to configure the pipeline for one's own running environment.
  
## Classifying ASVs:

`part2_illumina_dada2.pl` runs the DADA2 ASV classification pipeline on one or more runs in the same project. On a RHEL7 IGS machine:

`screen`  
`qlogin -P jravel-lab -l mem_free=500M -q interactive.q`  
 `export LD_LIBRARY_PATH=/usr/local/packages/gcc/lib64`  
`source /usr/local/packages/usepackage/share/usepackage/use.bsh`  
`use python-2.7`  
`use qiime`  
`cd <project directory created by part1.sh>`  
`part2_illumina_dada2.pl -i <comma-separated-input-run-names> -v <variable-region> -p <project-ID>`  

You can close the terminal at any time. To return to this process in another session, run:
`screen -r`

# To use the lastest version of this pipeline:
## Terminal
Choose one of the following options.

-To create a local git clone:

`git clone https://github.com/jbholm/IGS_dada2_pipeline.git`

-To update your local git clone while discarding your local development history:

`git fetch origin master`  
`git reset --hard FETCH_HEAD`

-To merge the changes into your local git clone:

`git checkout master`  
`git pull -X patience # Git may not be able to merge all the changes without manual editing`

# To use an earlier version of this pipeline:
## Desktop
Find the desired version on https://github.com/jbholm/IGS_dada2_pipeline/commits/master. Click the "**<>**" button to browse the repository. Click the green "**Clone or download**" dropdown and click "**Download ZIP**" to download this version of the pipeline.
## Terminal
In a local clone:  
Enter the following commands to save any existing changes to Git's stash and checkout a previous commit (change **1.0.0** to your desired version).  
`git add .`  
`git stash # To discard changes instead of stashing them, git reset --hard`  
`git rev-list master --grep='v1.0.0' | head -n 1 | xargs git checkout`  

To create a local clone:  
Follow the desktop instructions above to find the desired version of the pipeline. Click the green "**Clone or download**" dropdown and copy the HTTPS URL. Then run:
`git clone PASTE/HTTPS/URL/HERE`

# License

The contents of `taxonomy/silva/` are derived from the Silva database, and reformatted for DADA2 by Benjamin Callahan. They are available under the Silva dual-licensing model for academic and commercial users: https://www.arb-silva.de/silva-license-information/  
  
This repository's copy of the UNITE database at `taxonomy/sh_general_release_dynamic_01.12.2017.fasta` is available under [Attribution-ShareAlike (CC BY-SA)](https://creativecommons.org/licenses/by-sa/4.0/).
  
All other files in this respository are available under [GPLv3](https://github.com/jbholm/IGS_dada2_pipeline/blob/master/LICENSE).

