#!/bin/sh

use () 
{ 
    eval `/usr/local/packages/usepackage-1.13/bin/usepackage -b $*`
}

DEBUG=""
DBG=""
DRY_RUN=""
SKIP_ERR_THLD=""
FOR=""
REV=""
MAXN=''
MAXEE=""
TRUNCQ=""
RMPHIX=""
MAXLEN=""
MINQ=""
ONESTEP=""
PARAMS=""
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do
  case "$1" in
    -i) # update the Perl and bash documentation!!!
        RAW_PATH=$2
        shift 2
        ;;
    -r1|--r1)
        R1=$2
        shift 2
        ;;
    -r2|--r2)
        R2=$2
        shift 2
        ;;
    -i1|--i1)
        I1=$2
        shift 2
        ;;
    -i2|--i2)
        I2=$2
        shift 2
        ;;
    -p)
        PROJECT=$2
        shift 2
        ;;
    -r)
        RUN_ID=$2
        shift 2
        ;;
    -m)
        MAP=$2
        shift 2
        ;;
    -v)
        VAR=$2
        shift 2
        ;;
    --help|-h)
        perldoc "${0%/*}/pipeline.sh"
        exit    
        ;;
    --debug)
        DEBUG=$1
        shift 1
        ;;
    -dbg)
        DBG=$1
        shift 1
        ;;
    --dry-run)
        DRY_RUN=$1
        shift 1
        ;;
    --skip-err-thld)
        SKIP_ERR_THLD=$1
        shift 1
        ;;
    --dada2-truncLen-f|--for)
        FOR="$1 $2"
        shift 2
        ;;
    --dada2-truncLen-R|--rev)
        REV="$1 $2"
        shift 2
        ;;
    --dada2-maxN)
        MAXN="$1 $2"
        shift 2
        ;;
    --dada2-maxEE)
        MAXEE="$1 $2"
        shift 2
        ;;
    --dada2-truncQ)
        TRUNCQ="$1 $2"
        shift 2
        ;;
    --dada2-rmPhix)
        RMPHIX=$1
        shift 1
        ;;
    --dada2-maxLen)
        MAXLEN="$1 $2"
        shift 2
        ;;
    --dada2-minLen)
        MINLEN="$1 $2"
        shift 2
        ;;
    --dada2-minQ)
        MINQ="$1 $2"
        shift 2
        ;;
    --1Step)
        ONESTEP=$1
        shift 1
        ;;
    --storage-dir|--sd|-sd)
        SD=$2
        shift 2
        ;;
    -qp|--qp)
        QP=$2
        shift 2
        ;;
    *) # preserve positional arguments even if they fall between other params
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done
# Normally PARAMS would contain positional parameters, but our script doesn't take any
if [[ -n $PARAMS ]]; then
    echo "Passing unknown parameters through to illumina_dada2.pl: $PARAMS"
fi

# Acquire binaries
use sge
cd /usr/local/packages/qiime-1.9.1 || exit $?
. ./activate.sh
export PATH=/usr/local/packages/python-2.7.14/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/packages/python-2.7.14/lib:/usr/local/packages/gcc/lib64:$LD_LIBRARY_PATH
. /usr/local/packages/usepackage/share/usepackage/use.bsh

# Housekeeping
if [ "$SD" = "scratch" ]; then
    SD="/local/scratch/"
elif [ "$SD" = "groupshare" ]; then
    SD="/local/groupshare/ravel"
fi
DIR="$SD/$PROJECT/$RUN_ID/"
mkdir -p "$DIR/qsub_error_logs/"
mkdir -p "$DIR/qsub_stdout_logs/"

CMD=("-cwd" "-b y" "-l mem_free=200M" "-P jravel-lab" "-q threaded.q" "-pe thread 4" "-V" "-o ${DIR}qsub_stdout_logs/illumina_dada2.pl.stdout" "-e ${DIR}qsub_error_logs/illumina_dada2.pl.stderr" "/home/jolim/IGS_dada2_pipeline/illumina_dada2.pl" "-r1 $R1" "-r2 $R2" "-r3 $R3" "-r4 $R4" "-p $PROJECT" "-r $RUN_ID" "-sd $SD" "-v $VAR" "-m $MAP" "$DEBUG $DBG" "$DRY_RUN" "$SKIP_ERR_THLD" "$FOR" "$REV" "$MAXN" "$MAXEE" "$TRUNCQ" "$RMPHIX" "$MAXLEN" "$MINLEN" "$MINQ" "$ONESTEP" "$PARAMS")

echo "qsub ${CMD[*]}"
qsub ${CMD[*]}

: <<=cut
=pod

=head1 NAME

part1.sh

=head1 DESCRIPTION

The script can be launched from any location on the IGS server, it automatically 
produces a directory in /local/groupshare/ravel named after the project and run
ID provided. 

Beginning with the path to four raw Illumina sequencing files (R1, R2, I1, I2),
a mapping file, a project ID, a run ID, and specifying the targeted variable
region, this script:
  1. Produces individual .fastq files for each sample listed in the mapping file
  2. Performs tag-cleaning of each file
  3. Runs the R1 (forward) and R2 (reverse) files through the dada2 pipeline for either the 16S rRNA gene V3V4 or V4 regions.

A log file is written at <PROJECT>/<RUN>/<PROJECT>_<RUN>_16S_pipeline_log.txt

=head1 SYNOPSIS

THE FOLLOWING STEPS ARE REQUIRED BEFORE THIS SCRIPT WILL RUN:

1. Qlogin to RHEL7 if you ssh'd onto a non-RHEL7 host: 

  qlogin -P jravel-lab -l mem_free=500M -q interactive.q -V

2. Execute from bash: 

  export LD_LIBRARY_PATH=/usr/local/packages/gcc/lib64
  source /usr/local/packages/usepackage/share/usepackage/use.bsh
  use python-2.7
  use qiime

3. then run script as below: 

FOR 2-STEP:
part1.sh -i <path to raw files> -p <project name> -r <run id> 
  -v <variable region> -m <full path to mapping file> -sd <storage directory>

OR:
part1.sh -r1 <path_to_R1_file> -r2 <path_to_R2_file> -i1 <path_to_I1_file>
  -i2 <path_to_I2_file> -p <project name> -r <run id> -v <variable region> 
  -m <full path to mapping file> -sd <storage directory>

FOR 1-STEP:
part1.sh -i <input directory> -p <project name> -r <run ID> -m <mapping file>
  -v <variable region> -sd <storage directory> --1Step

OR:
part1.sh -r1 <path_to_R1_file> -r2 <path_to_R2_file> -p <project name> 
  -r <run ID> -m <mapping file> -v <variable region> -sd <storage directory> 
  --1Step

=head1 OPTIONS

=over

=item B<--raw-path>=path, B<-i> path

Single full path to directory containing raw R1, R2, R3, and R4 files

=item B<--r1-path>=file, B<-r1> file

Full path to raw R1 read file (forward read file, or r1) (.fastq.gz).

=item B<--r2-path>=file, B<-r2> file

Full path to raw R2 read file (barcode file, or i1) (.fastq.gz).

=item B<--r3-path>=file, B<-r3> file

Full path to raw R3 read file (barcode file, or i2) (.fastq.gz).

=item B<--r4-path>=file, B<-r4> file

Full path to raw R4 read file (reverse read file, r2 or r4) (.fastq.gz).

=item B<--project-name>=name, B<-p> name

Create the project folder with this name.

=item B<--run-ID>=name, B<-r> name

Create the run folder with this name.

=item B<--map>=file, B<-m> file

The full path to the Qiime-formatted mapping file.

=item B<--var-reg>={V3V4, V4, ITS}, B<-v> {V3V4, V4, ITS}

The targeted variable region.

=item B<--1Step>

Use this flag if the data are prepared by 1-Step PCR (only r1 & r2 raw files
available)

=item B<--storage-dir>=path, B<-sd> path

Indicate an existing directory in which to place the project directory.
"scratch" evaluates to "/local/scratch/" and "groupshare" evaluates to
"/local/groupshare/ravel"

=item B<-h>, B<--help>

Print help message and exit successfully.

=item B<--qsub-project>=space, B<-qp> space

Indicate which qsub-project space should be used for all qsubmissions. The
default is jravel-lab.

=item B<-dbg> {qiime_and_validation, extract_barcodes, demultiplex, tagclean, dada2}

Runs only one section of the pipeline and writes every qsub command to the log file.

=back 

=cut

