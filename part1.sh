#!/bin/bash

use () 
{ 
    eval `/usr/local/packages/usepackage-1.13/bin/usepackage -b "$*"`
}

QSUB_ARGS=""
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

while [[ ! "$1" == "--" && "$#" != 0 ]]; do
  case "$1" in
    --qsub)
        QSUB_ARGS="${1#*=}"
        shift 1
        ;;
    -i) # update the Perl and bash documentation!!!
        RAW_PATH="-i $2"
        shift 2
        ;;
    -r1|--r1)
        R1="$2"
        shift 2
        ;;
    -r2|--r2)
        R2="$2"
        shift 2
        ;;
    -i1|--i1)
        I1="$2"
        shift 2
        ;;
    -i2|--i2)
        I2="$2"
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
        perldoc -F "${0}"
        exit 0
        ;;
    -dbg)
        DBG="$DBG -dbg $2"
        shift 2
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
    # note that options and their operands will be separate elements of $PARAMS
      PARAMS="$PARAMS $1"
      shift 1
      ;;
  esac
done

# Normally PARAMS would contain positional parameters, but our script doesn't take any
if [[ -n $PARAMS ]]; then
    printf "\nPassing unknown parameters through to illumina_dada2.pl: $PARAMS\n"
fi

# validate -dbg flags
if [[ -n "$DBG" ]]; then # if dbg is non-empty string
    DBG=${DBG:1:${#DBG}} # Remove the first character (whitespace) from dbg
    DBG_A=( $DBG ) # Separate on whitespace and convert to array

    for ((WORD=7;WORD<${#DBG_A[@]};WORD++)); do
        if (( WORD % 2 == 0 )); then
            if [ ${DBG_A[WORD]} != "validate" ] && \
            [ ${DBG_A[WORD]} != "barcodes" ] && \
            [ ${DBG_A[WORD]} != "demultiplex" ] && \
            [ ${DBG_A[WORD]} != "tagclean" ] && \
            [ ${DBG_A[WORD]} != "dada2" ]; then
                printf "\nIllegal debug option. Legal debug options are validate, "\
                "barcodes, demultiplex, tagclean, and dada2.\n"
                exit 2
            fi
        fi
    done
    NDBG="${#DBG_FLAGS[@]}" # Count the number of dbg flags
    printf "%b" "Debug options validated\n\n"
fi

# Validate specification of raw files / directory
if [[ ! -n "$RAW_PATH" ]]; then
    if [[ -n "$ONESTEP" ]]; then
        if ! [[ -n "$R1" && -n "$R2" ]]; then
            printf "%b\n" "Please provide the location of the raw sequencing files (single directory => -i)" \
            "\tOR" \
            "Full paths to each raw file. A 1-step run requires -r1 and -r2.\n"
            exit 2
        elif ! [[ -e "$R1" && -e "$R2" ]]; then
            printf "%b\n" "Unable to find -r1, -r2. Please check spellings and " \
            "file permissions.\n"
            exit 2
        fi
        R1="-r1 $R1"
        R2="-r2 $R2"
    else
        if ! [[ -n "$R1" && -n "$R2" && -n "$I1" && -n "$I2" ]]; then
            printf "%b\n" "Please provide the location of the raw sequencing files (single directory => -i)" \
            "\tOR" \
            "Full paths to each raw file. A 2-step run requires -r1, -r2, -i1, and -i2.\n"
            exit 2
        elif ! [[ -e "$R1" && -e "$R2" && -e "$I1" && -e "$I2" ]]; then
            printf "%b\n" "Unable to find -r1, -r2, -i1, and -i2. Please check " \
            "spellings and file permissions.\n"
            exit 2
        fi
        R1="-r1 $R1"
        R2="-r2 $R2"
        I1="-i1 $I1"
        I2="-i2 $I2"
    fi
elif [[ -n "$R1" || -n "$R2" || -n "$I1" || -n "$I2" ]]; then
    printf "%b\n" "Do not provide both the input directory (-i) and the input files." \
    "(-r1, -r2, -i1, -i2). Only one or the other is needed.\n"
    exit 2
elif [[ ! -d "$RAW_PATH" ]]; then
    printf "\n\tThe input directory does not exist or is inaccessible.\n"
    exit 2
fi

# Validate that mapping file was given and exists
if [[ ! -n "$MAP" ]]; then
    printf "\n\tPlease provide a full path to the project mapping file (-m)\n\n"
    exit 2
elif [[ ! -e "$MAP" ]]; then
    printf "\n\tMapping file (-m) does not exist or is inaccessible.\n\n"
    exit 2
fi

# Validate the variable region
if [[ ! -n "$VAR" ]]; then
    printf "%b" "Please indicate the targeted variable region (-v V3V4 or -v V4" \
    " or -v ITS)\n\n"
    exit 2
elif [[ !( "$VAR" == "V3V4" || "$VAR" == "V4" || "$VAR" == "ITS" ) ]]; then
    printf "\n\tVariable region was '$VAR' but only 'V3V4', 'V4', and 'ITS' are"\
    " supported.\n\n"
    exit 2
fi

if [[ ! -n "$SD" ]]; then
    printf "%b" "***Please choose a storage directory (-sd), either 'scratch' or "\
    "'groupshare'.\n\tOR\nprovide a full path to an existing directory.\n\n"
    exit 2
fi

if [[ ( -n $FOR && ! -n $REV ) || ( -n $REV && ! -n $FOR ) ]]; then
    printf "***\nPlease provide truncation lengths for forward and reverse "\
    "reads\n";
    exit 2
fi

if [[ -n "$qproj" ]]; then
    printf "\nqsub-project ID (--qsub-project, -qp) not provided. Using "\
    "jravel-lab as default\n"
    $qproj = "jravel-lab";
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
    SD="/local/groupshare/ravel/"
fi
DIR="${SD}${PROJECT}/${RUN_ID}/"
mkdir -p "$DIR/qsub_error_logs/"
mkdir -p "$DIR/qsub_stdout_logs/"

# Begin log (will be continued by illumina_dada2.pl)
# Print pipeline version, system time, and qsub/illumina_dada2.pl command
log="$SD/${PROJECT}_${RUN_ID}_16S_pipeline_log.txt";
MY_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
LOG_VERSION="$MY_DIR/scripts/log_version.pl"
perl $LOG_VERSION $log
printf "$time\n" > $log 

CMD=("$QSUB_ARGS" "-cwd" "-b y" "-l mem_free=200M" "-P jravel-lab" "-q threaded.q" "-pe thread 4" "-V" "-o ${DIR}qsub_stdout_logs/illumina_dada2.pl.stdout" "-e ${DIR}qsub_error_logs/illumina_dada2.pl.stderr" "/home/jolim/IGS_dada2_pipeline/illumina_dada2.pl" "$RAW_PATH" "$R1" "$R2" "$I1" "$I2" "-p $PROJECT" "-r $RUN_ID" "-sd $SD" "-v $VAR" "-m $MAP" "$DEBUG $DBG" "$DRY_RUN" "$SKIP_ERR_THLD" "$FOR" "$REV" "$MAXN" "$MAXEE" "$TRUNCQ" "$RMPHIX" "$MAXLEN" "$MINLEN" "$MINQ" "$ONESTEP" "$PARAMS")

printf "Initial qsub command: \n$ qsub ${CMD[*]}\n" > $log
#qsub ${CMD[*]}

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

=item B<--qsub>="options"

Adds options to the outermost qsub call. By default, qsub is called with the
following options:

    -cwd
    -b y
    -l mem_free=200M 
    -P jravel-lab 
    -q threaded.q 
    -pe thread 4 
    -V
    -o <path auto-generated from -sd, -p, and -r>
    -e <from auto-generated path -sd, -p, and -r>

Additional options may be specified as a single string surrounded by
double quotes ("), as shown below:

    part1.sh --qsub="-m ea -l excl=true" -p project -r run -sd scratch...

=back 

=cut

