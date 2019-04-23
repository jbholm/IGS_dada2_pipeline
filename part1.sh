#!/bin/bash
set -e # exit when any command fails (returns non-0 exit status)

use () 
{ 
    eval `/usr/local/packages/usepackage-1.13/bin/usepackage -b "$*"`
}

# QSUB_ARGS=""
# DRY_RUN=""
# FOR=""
# REV=""
# MAXN=''
# MAXEE=""
# TRUNCQ=""
# RMPHIX=""
# MAXLEN=""
# MINQ=""
# ONESTEP=""
# PARAMS=""

try_assign()
{
    if [[ "$3" =~ ^- || ! -n "$3" ]]; then
        printf "%b" "$2 missing its value. Unable to continue.\n\n"
        exit 2
    else 
        eval "$1='$3'"
    fi
    return 0
}

while [[ ! "$1" == "--" && "$#" != 0 ]]; do
  case "$1" in
    --qsub=*)
        QSUB_ARGS="${1#*=}"
        if [[ ! -n "$QSUB_ARGS" ]]; then  
            printf "%b" "--qsub missing value. --qsub=\"\" and --qsub= are " \
            "not accepted.\n"
            exit 2
        fi
        shift 1
        ;;
    -i) # update the Perl and bash documentation!!!
        try_assign RAW_PATH "$1" "$2"
        shift 2
        ;;
    -r1|--r1)
        try_assign R1 "$1" "$2"
        shift 2
        ;;
    -r2|--r2)
        try_assign R2 "$1" "$2"
        shift 2
        ;;
    -i1|--i1)
        try_assign I1 "$1" "$2"
        shift 2
        ;;
    -i2|--i2)
        try_assign I2 "$1" "$2"
        shift 2
        ;;
    -p)
        try_assign PROJECT "$1" "$2"
        shift 2
        ;;
    -r)
        try_assign RUN "$1" "$2"
        shift 2
        ;;
    -m)
        try_assign MAP "$1" "$2"
        shift 2
        ;;
    -v)
        try_assign VAR "$1" "$2"
        shift 2
        ;;
    --help|-h)
        perldoc -F "${0}"
        exit 0
        ;;
    -dbg)
        if [[ "$2" =~ ^- || ! -n "$2" ]]; then
            printf "%b" "-dbg missing its value. Unable to continue.\n\n"
            exit 2
        else 
            DBG="$DBG -dbg $2"
            shift 2
        fi
        ;;
    --dry-run)
        DRY_RUN=$1
        shift 1
        ;;
    --verbose)
        VERBOSE=$1
        shift 1
        ;;
    --skip-err-thld)
        SKIP_ERR_THLD=$1
        shift 1
        ;;
    --dada2-truncLen-F|--for)
        try_assign FOR "$1" "$2"
        shift 2
        ;;
    --dada2-truncLen-R|--rev)
        try_assign REV "$1" "$2"
        shift 2
        ;;
    --dada2-maxN)
        try_assign MAXN "$1" "$2"
        shift 2
        ;;
    --dada2-maxEE)
        try_assign MAXEE "$1" "$2"
        shift 2
        ;;
    --dada2-truncQ)
        try_assign TRUNCQ "$1" "$2"
        shift 2
        ;;
    --dada2-rmPhix)
        try_assign RMPHIX "$1" "$2"
        shift 2
        ;;
    --dada2-maxLen)
        try_assign MAXLEN "$1" "$2"
        shift 2
        ;;
    --dada2-minLen)
        try_assign MINLEN "$1" "$2"
        shift 2
        ;;
    --dada2-minQ)
        try_assign MINQ "$1" "$2"
        shift 2
        ;;
    --1Step)
        ONESTEP=$1
        shift 1
        ;;
    --storage-dir|-sd)
        try_assign SD "$1" "$2"
        shift 2
        ;;
    -qp|--qp)
        try_assign QP "$1" "$2"
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

if [[ -n "$ONESTEP" ]]; then
    printf "RUN TYPE: One-step\n"
else
    printf "RUN TYPE: Two-step\n"
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
        printf "FWD READS: $R1\nREV READS:$R2\n"
        INPUT="-r1 $R1 -r2 $R2"
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
        printf "%b\n" "FWD READS: $R1" \
            "REV READS: $R2" \
            "INDEX 1: $I1" \
            "INDEX 2: $I2"
        INPUT="-r1 $R1 -r2 $R2 -i1 $I1 -i2 $I2"
    fi
else
    if [[ -n "$R1" || -n "$R2" || -n "$I1" || -n "$I2" ]]; then
        printf "%b\n" "Do not provide both the input directory (-i) and the input files." \
        "(-r1, -r2, -i1, -i2). Only one or the other is needed.\n"
        exit 2
    elif [[ ! -d "$RAW_PATH" ]]; then
        printf "\n\tThe input directory does not exist or is inaccessible.\n"
        exit 2
    else
        printf "INPUT DIRECTORY: $RAW_PATH\n"
    fi
    INPUT="-i $RAW_PATH"
fi

# Validate that mapping file was given and exists
if [[ ! -n "$MAP" ]]; then
    printf "\n\tPlease provide a full path to the project mapping file (-m)\n\n"
    exit 2
elif [[ ! -e "$MAP" ]]; then
    printf "\n\tMapping file (-m) does not exist or is inaccessible.\n\n"
    exit 2
fi
printf "MAPPING FILE: $MAP\n"

# -p is mandatory
if [[ ! -n "$PROJECT" ]]; then
    printf "-p <project> required.\n\n"
    exit 2
fi

# -r is mandatory
if [[ ! -n "$RUN" ]]; then
    printf "-r <run ID> required.\n\n"
    exit 2
fi

# Validate the variable region
if [[ ! -n "$VAR" ]]; then
    VAR="V3V4"
elif [[ !( "$VAR" == "V3V4" || "$VAR" == "V4" || "$VAR" == "ITS" ) ]]; then
    printf "%b" "\nVariable region was '$VAR' but only 'V3V4', 'V4', and 'ITS'"\
    " are supported.\n\n"
    exit 2
fi
printf "VARIABLE REGION: $VAR\n"

# Validate the storage directory
if [[ ! -n "$SD" || "$SD" == "scratch" ]]; then
    SD="/local/scratch"
elif [[ "$SD" == "groupshare" ]]; then
    SD="/local/groupshare/ravel"
elif [[ ! -d "$SD" ]]; then
    printf "$SD does not exist!\n" # A custom storage directory must exist
    exit 2
else
    SD="${SD%\/}" # Remove any trailing slash
fi
SD="$SD/$PROJECT/$RUN" # But we'll automatically make these subfolders
printf "%b" "WORKING DIRECTORY: $SD\n"
mkdir -p "$SD/qsub_error_logs/"
mkdir -p "$SD/qsub_stdout_logs/"

# validate -dbg flags
if [[ -n "$DBG" ]]; then # if dbg is non-empty string
    DBG=${DBG:1:${#DBG}} # Remove the first character (whitespace) from dbg
    DBG_A=( $DBG ) # Separate on whitespace and convert to array

    for ((WORD=7;WORD<${#DBG_A[@]};WORD++)); do
        if (( WORD % 2 == 0 )); then
            if [ ${DBG_A[WORD]} != "validate" ] && \
            [ ${DBG_A[WORD]} != "barcodes" ] && \
            [ ${DBG_A[WORD]} != "demux" ] && \
            [ ${DBG_A[WORD]} != "tagclean" ] && \
            [ ${DBG_A[WORD]} != "dada2" ]; then
                printf "\nIllegal debug option. Legal debug options are validate, "\
                "barcodes, demux, tagclean, and dada2.\n"
                exit 2
            else
                printf "DEBUG: ${DBG_A[WORD]}\n"
            fi
        fi
    done
    #NDBG="${#DBG_FLAGS[@]}" # Count the number of dbg flags
fi

if [[ -n "$VERBOSE" ]]; then
    printf "%b\n" "Running verbose..." \
    "All qsub commands will be printed to: ${SD}/qsub_stdout_logs/illumina_dada2.pl.stdout"
fi

if [[ ( -n $FOR && ! -n $REV ) || ( -n $REV && ! -n $FOR ) ]]; then
    printf "***\nPlease provide truncation lengths for forward and reverse "\
    "reads\n";
    exit 2
fi

# Acquire binaries
use sge
cd /usr/local/packages/qiime-1.9.1 || exit $?
. ./activate.sh
export PATH=/usr/local/packages/python-2.7.14/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/packages/python-2.7.14/lib:/usr/local/packages/gcc/lib64:$LD_LIBRARY_PATH
. /usr/local/packages/usepackage/share/usepackage/use.bsh

# Begin log (will be continued by illumina_dada2.pl)
# Print pipeline version, system time, and qsub/illumina_dada2.pl command
log="$SD/${PROJECT}_${RUN}_16S_pipeline_log.txt";
MY_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
LOG_VERSION="$MY_DIR/scripts/log_version.pl"
perl $LOG_VERSION $log
printf "`date`\n" >> $log 

CMD=("$QSUB_ARGS" "-cwd" "-b y" "-l mem_free=200M" "-P jravel-lab" "-q threaded.q" "-pe thread 4" "-V" "-o ${SD}/qsub_stdout_logs/illumina_dada2.pl.stdout" "-e ${SD}/qsub_error_logs/illumina_dada2.pl.stderr" "/home/jolim/IGS_dada2_pipeline/illumina_dada2.pl" "$INPUT" "-wd $SD" "-v $VAR" "-m $MAP" "$DBG" "$VERBOSE" "$DRY_RUN" "$SKIP_ERR_THLD" "$FOR" "$REV" "$MAXN" "$MAXEE" "$TRUNCQ" "$RMPHIX" "$MAXLEN" "$MINLEN" "$MINQ" "$ONESTEP" "$PARAMS")
printf "$ qsub${CMD[*]}\n"
if [[ -n $VERBOSE ]]; then
    printf "$ qsub${CMD[*]}\n" >> $log
fi
qsub ${CMD[*]}

: <<=cut
=pod

=head1 NAME

part1.sh

=head1 DESCRIPTION

The script can be launched from any location on the IGS server, it automatically 
produces a directory in /local/groupshare/ravel named after the project and run
ID provided. 

Beginning with a set of raw Illumina sequencing files (usually R1, R2, I1, and I2),
a mapping file, a project ID, a run ID, and specifying the targeted variable
region, this script:
  1. Extracts barcodes from the raw files.
  2. Demultiplexes the raw reads into fastq files containing reads specific to 
  this project.
  3. Produces individual .fastq files for each sample listed in the mapping file
  4. Performs tag-cleaning of each file
  5. Runs the forward and reverse reads through the dada2 pipeline for the 
  specified 16S rRNA gene region.

A log file is written at <PROJECT>/<RUN>/<PROJECT>_<RUN>_16S_pipeline_log.txt


=head1 SYNOPSIS

part1.sh [-i <input directory> | -r1 <fwd reads> -r2 <rev reads> [-i1 <index 1> -i2 <index 2>]] -p <project> -r <run> -m <map> [-v <variable region>] [--1Step] [<options>]

=head1 OPTIONS

=over

=item B<--raw-path>=path, B<-i> path

Single full path to directory containing raw files. File names must have the 
pattern *AA.fastq[.gz], where AA is either R1, R2, R3, R4, I1, or I2, and gzip 
compression is optional. B<--raw-path> is incompatible with B<-r1>, B<-r2>, B<-i1>, and B<-i2>.

Three combinations of files are allowed:

    DEFAULT (TWO-STEP PCR):
    R1, R2, I1, I2

    OLD TWO-STEP NAMING STYLE (automatically detected):
    R1, R2, R3, R4 
    (R1 and R4 are fwd and rev reads, respectively. R2 and R3 are index 1 
    and index 2, respectively.)

    ONE-STEP PCR (<--1step> MUST BE GIVEN):
    R1, R2

=item B<-r1> file

Full path to raw foward read file (R1). Gzip compression optional.

=item B<-r2> file

Full path to raw reverse read file (R2, or R4 in old naming scheme). Gzip 
compression optional.

=item B<-i1> file

Full path to raw index 1 file (I1, or R2 in old naming scheme). Incompatible 
with B<--1step>. Gzip compression optional.

=item B<-i2> file

Full path to raw index 2 file (I2, or R3 in old naming scheme). Incompatible 
with B<--1step>. Gzip compression optional.

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

=item B<-dbg> {validate, barcodes, demux, tagclean, dada2}

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
double quotes ("), as shown below. Most options will override this script's defaults.

    part1.sh --qsub="-m ea -l excl=true" -p project -r run -sd /other/path...

=back 

=cut

