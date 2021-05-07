#!/bin/bash
set -e # exit when any command fails (returns non-0 exit status)

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

stop()
{
    printf "\n%b\n\n" "$1"
    exit 2
}

try_assign()
{
    if [[ "$3" =~ ^- || ! -n "$3" ]]; then
        stop "$2 missing its value. Unable to continue."
    else 
        eval "$1='$3'"
    fi
    return 0
}

assign_default()
{
    if [[ -n "$1" ]]; then
        eval "$1='$2'"
        return 1
    fi
    return 0
}

while [[ ! "$1" == "--" && "$#" != 0 ]]; do
  case "$1" in
    --qsub*)
        QSUB_ARGS="${1#*=}"
        if [[ ! -n "$QSUB_ARGS" || ! $1 =~ "=" ]]; then  
            MSG="--qsub missing value."
            MSG+=" --qsub=\"\" and --qsub= are not accepted."
            stop "$MSG"
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
    -p|--project-name)
        try_assign PROJECT "$1" "$2"
        shift 2
        ;;
    -r|--run-ID)
        try_assign RUN "$1" "$2"
        shift 2
        ;;
    -m|--map)
        try_assign MAP "$1" "$2"
        shift 2
        ;;
    --bclen)
        if [[ "$2" =~ ^- || ! -n "$2" ]]; then
            stop "--bclen missing its value. Unable to continue."
        else 
            BCLENGTH="--bclen $2"
            shift 2
        fi
        ;;
    --troubleshoot_barcodes)
        TROUBLESHOOT_BARCODES=$1
        shift 1
        ;;
    -v|--var-reg)
        try_assign VAR "$1" "$2"
        shift 2
        ;;
    --help|-h)
        perldoc -F "${0}"
        exit 0
        ;;
    --debug|-d)
        if [[ "$2" =~ ^- || ! -n "$2" ]]; then
            stop "--debug missing its value. Unable to continue."
        else 
            DBG="$DBG --debug $2"
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
    --dada2-mem)
        if [[ "$2" =~ ^- || ! -n "$2" ]]; then
            stop "--debug missing its value. Unable to continue."
        else 
            DADA2MEM="--dada2-mem $2"
            shift 2
        fi
        ;;
    --dada2*)
        DADA2="$1"
        if [[ ! -n "${DADA2#*--dada2}" || ! $DADA2 =~ "=" ]]; then  
            MSG="--dada2 missing value."
            MSG+=" --dada2=\"\" and --dada2= are not accepted."
            stop "$MSG"
        fi
        shift 1
        ;;
    --1Step)
        ONESTEP=$1
        shift 1
        ;;
    --storage-dir|-sd)
        try_assign SD "$1" "$2"
        shift 2
        ;;
    -qp|--qsub-project)
        try_assign QP "$1" "$2"
        shift 2
        ;;
    --email)
        EMAIL="-m ea"
        shift 1
        ;;
    *) # preserve positional arguments even if they fall between other params
    # note that options and their operands will be separate elements of $PARAMS
      PARAMS="$PARAMS $1"
      shift 1
      ;;
  esac
done

# directory of this script (part1.sh)
MY_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

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
            MSG="Please provide the location of the raw sequencing files"
            MSG+=" (single directory => -i)\n"
            MSG+="\tOR\n"
            MSG+="Full paths to each raw file. "
            MSG+="A 1-step run requires -r1 and -r2."
            stop "$MSG"
        elif ! [[ -e "$R1" && -e "$R2" ]]; then
            MSG="Unable to find -r1, -r2. Please check spellings and "
            MSG+="file permissions."
            stop "$MSG"
        fi
        printf "FWD READS: $R1\nREV READS:$R2\n"
        INPUT="-r1 $R1 -r2 $R2"
    else
        if ! [[ -n "$R1" && -n "$R2" && -n "$I1" && -n "$I2" ]]; then
            MSG="Please provide the location of the raw sequencing files "
            MSG+="(single directory => -i)\n"
            MSG+="\tOR\n"
            MSG+="Full paths to each raw file. "
            MSG+="A 2-step run requires -r1, -r2, -i1, and -i2."
            stop "$MSG"
        elif ! [[ -e "$R1" && -e "$R2" && -e "$I1" && -e "$I2" ]]; then
            MSG="Unable to find -r1, -r2, -i1, and -i2. Please check "
            MSG+="spellings and file permissions."
            stop "$MSG"
        fi
        printf "%b\n" "FWD READS: $R1" \
            "REV READS: $R2" \
            "INDEX 1: $I1" \
            "INDEX 2: $I2"
        INPUT="-r1 $R1 -r2 $R2 -i1 $I1 -i2 $I2"
    fi
else
    if [[ -n "$R1" || -n "$R2" || -n "$I1" || -n "$I2" ]]; then
        MSG="Do not provide both the input directory (-i) and the input "
        MSG+="files (-r1, -r2, -i1,\n-i2). Only one or the other is needed."
        stop "$MSG"
    elif [[ ! -d "$RAW_PATH" ]]; then
        stop "\tThe input directory does not exist or is inaccessible."
    else
        printf "INPUT DIRECTORY: $RAW_PATH\n"
    fi
    INPUT="-i $RAW_PATH"
fi

# Validate that mapping file was given and exists
if [[ ! -n "$MAP" ]]; then
    stop "\tPlease provide a full path to the project mapping file (-m)"
elif [[ ! -e "$MAP" ]]; then
    stop "\tMapping file (-m) does not exist or is inaccessible."
fi
printf "MAPPING FILE: $MAP\n"

# -p is mandatory
if [[ ! -n "$PROJECT" ]]; then
    stop "Project name (-p) required."
fi

# -r is mandatory
if [[ ! -n "$RUN" ]]; then
    stop "Run ID (-r) required."
fi

# Validate the variable region
if [[ ! -n "$VAR" ]]; then
    VAR="V3V4"
elif [[ !( "$VAR" == "V3V4" || "$VAR" == "V4" || "$VAR" == "ITS" ) ]]; then
    MSG="Variable region was '$VAR' but only 'V3V4', 'V4', and 'ITS' are "
    MSG+="supported."
    stop "$MSG"
fi
printf "VARIABLE REGION: $VAR\n"

# Validate the storage directory
if [[ ! -n "$SD" || "$SD" == "scratch" ]]; then
    SD="/local/scratch"
elif [[ "$SD" == "groupshare" ]]; then
    SD="/local/groupshare/ravel"
elif [[ ! -d "$SD" ]]; then
    stop "$SD does not exist!\n" # A custom storage directory must exist
else
    SD="${SD%\/}" # Remove any trailing slash
fi
SD="$SD/$PROJECT/$RUN" # But we'll automatically make these subfolders
printf "%b" "WORKING DIRECTORY: $SD\n"
mkdir -p "$SD/qsub_error_logs/"
mkdir -p "$SD/qsub_stdout_logs/"

if [[ ! -n "$QP" ]]; then
    printf "qsub-project ID (--qp) not provided. Using jravel-lab as default\n"
    QP="jravel-lab"
fi

# validate -dbg flags
if [[ -n "$DBG" ]]; then # if dbg is non-empty string
    DBG=${DBG:1:${#DBG}} # Remove the first character (whitespace) from dbg
    DBG_A=( $DBG ) # Separate on whitespace and convert to array

    # It's allowed for the --debug flags to be out of order.
    # But, starting with the earliest step specified, only the first "streak"
    # of consecutive steps will be run.
    DBG=""
    CONSECUTIVE=false
    EARLIEST=true
    for STEP in "barcodes" "demux" "splitsamples" "tagclean" "dada2"; do
        PRESENT=false
        for ((WORD=0;WORD<${#DBG_A[@]};WORD++)); do
            if (( WORD % 2 == 1 )); then
                if [ ${DBG_A[WORD]} != "barcodes" ] && \
                [ ${DBG_A[WORD]} != "demux" ] && \
                [ ${DBG_A[WORD]} != "splitsamples" ] && \
                [ ${DBG_A[WORD]} != "tagclean" ] && \
                [ ${DBG_A[WORD]} != "dada2" ]; then
                    MSG="Illegal debug option ${DBG_A[WORD]}. Legal debug options are "
                    MSG+="barcodes, demux, splitsamples, tagclean, and dada2."
                    stop "$MSG"
                else
                    if [[ $STEP == ${DBG_A[WORD]} ]]; then
                        if [[ $EARLIEST = true || $CONSECUTIVE = true ]]; then
                            printf "DEBUG: $STEP\n"
                            # start the "streak" on the earliest flag found
                            CONSECUTIVE=true 
                            EARLIEST=false
                            # Pass along the debug flags only if they are part
                            # of the streak
                            DBG+="--debug $STEP "
                        fi
                        if [[ $EARLIEST = true ]]; then
                            EARLIEST=false
                        fi
                        PRESENT=true
                    fi
                fi
            fi
        done
        # If any step not present, break the "streak" of consecutive steps
        if [[ ! $PRESENT = true ]]; then
            CONSECUTIVE=false
        fi
    done

    printf "%b\n" "Use --debug flags with care!" \
    "Input files for your selected steps cannot be checked for consistency with the" \
    "raw Illumina files."
fi

if [[ -n "$VERBOSE" ]]; then
    printf "%b\n" "Running verbose..." \
    "All shell commands will be printed to: \n${SD}/qsub_stdout_logs/illumina_dada2.pl.stdout"
fi

if [[ -n "$EMAIL" ]]; then
    QSUB_ARGS="$QSUB_ARGS $EMAIL"
fi

# Acquire binaries

use () 
{ 
    eval `/usr/local/packages/usepackage/bin/usepackage -b $*` || true
}


use sge > /dev/null 2>&1  || true
module load sge 2>/dev/null || true

# if [[ -e "/usr/local/packages/miniconda3/etc/profile.d/conda.sh" ]]; then
#     source /usr/local/packages/miniconda3/etc/profile.d/conda.sh 2>/dev/null
# else
#     source "/local/projects-t3/MSL/pipelines/packages/miniconda3/etc/profile.d/conda.sh" 2>/dev/null
# fi
# conda=`cat "$MY_DIR/config.json" | \
#     python3 -c "import sys, json; print(json.load(sys.stdin)['conda'])"`
# qiime_env=`cat "$MY_DIR/config.json" | \
#     python3 -c "import sys, json; print(json.load(sys.stdin)['qiime_env'])"`
# conda activate "$qiime_env"

# . /usr/local/packages/qiime-1.9.1/activate.sh
# export PATH=/usr/local/packages/python-2.7/bin:$PATH
# export LD_LIBRARY_PATH=/usr/local/packages/python-2.7/lib:/usr/local/packages/gcc/lib64:$LD_LIBRARY_PATH
# . /usr/local/packages/usepackage/share/usepackage/use.bsh 2>/dev/null || true

# # Begin log (will be continued by illumina_dada2.pl)
log="$SD/${PROJECT}_${RUN}_16S_pipeline_log.txt"

# Remove extra spaces caused by joining empty arguments with a whitespace
OPTSARR=("$PARAMS" "$BCLENGTH" "$TROUBLESHOOT_BARCODES" "$ONESTEP" "$DADA2" "$DADA2MEM" "$DBG" "$VERBOSE" "$DRY_RUN")
OPTS="${OPTSARR[*]}"
OPTS="$( echo "$OPTS" | awk '{$1=$1;print}' )"

ARGS=("-cwd" "-b y" "-l mem_free=4G" "-P" "$QP" "-q threaded.q" "-pe thread 4" "-V" "-N" "MSL_$PROJECT" "-o ${SD}/qsub_stdout_logs/illumina_dada2.pl.stdout" "-e ${SD}/qsub_error_logs/illumina_dada2.pl.stderr" "$QSUB_ARGS" "${MY_DIR}/illumina_dada2.pl" "$INPUT" "-wd" "$SD" "-v" "$VAR" "-m" "$MAP" "$OPTS")
CMD=()
for ARG in "${ARGS[@]}"; do
    if [[ -n "$ARG" ]]; then
        CMD+=("$ARG")
    fi 
done

printf "$ qsub ${CMD[*]}\n"
printf "$ qsub ${CMD[*]}\n" >> $log

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
a mapping file, a project ID, and a run ID, this script:

1. Extracts barcodes from the raw files.

2. Demultiplexes the raw reads into fastq files containing reads specific to 
this project.

3. Produces individual .fastq files for each sample listed in the mapping file

4. Performs tag-cleaning of each sample-specific file

5. Runs the forward and reverse reads through the dada2 pipeline for the 
V3V4 16S rRNA gene region. Alternatively, analysis of the V4 or ITS region may
be specified.

A log file is written at: <PROJECT>/<RUN>/<PROJECT>_<RUN>_16S_pipeline_log.txt


=head1 SYNOPSIS

part1.sh (-i <input directory> | -r1 <fwd reads> -r2 <rev reads> [-i1 <index 1> -i2 <index 2>]) -p <project> -r <run> -m <map> [-v <variable region>] [--1Step] [<options>]

=head1 OPTIONS

=over

=item B<-i> path

Single full path to directory containing raw files. File names must have the 
pattern *AA.fastq[.gz], where * is wildcard, AA is either R1, R2, R3, R4, I1, or
I2, and gzip compression is optional. B<-i> is incompatible with B<-r1>, B<-r2>,
B<-i1>, and B<-i2>.

Three combinations of files are allowed:

    DEFAULT (TWO-STEP PCR):
    R1, R2, I1, I2

    OLD TWO-STEP NAMING STYLE (automatically detected):
    R1, R2, R3, R4 
    (R1 and R4 are fwd and rev reads, respectively. R2 and R3 are index 
    1 and index 2, respectively.)

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

=item B<--project-name> name, B<-p> name

Create the project folder with this name.

=item B<--run-ID> name, B<-r> name

Create the run folder with this name.

=item B<--map> file, B<-m> file

The full path to the Qiime-formatted mapping file.

=item B<--var-reg> {V3V4, V4, ITS}, B<-v> {V3V4, V4, ITS}

The targeted variable region. V3V4 is default.

=item B<--1Step>

Use this flag if the data are prepared by 1-Step PCR (only r1 & r2 raw files
available)

=item B<--storage-dir> path, B<-sd> path

Indicate an existing directory in which to place the project directory.
"scratch" (default) evaluates to "/local/scratch/" and "groupshare" evaluates to
"/local/groupshare/ravel".

=item B<--bclen> LENGTH

The length of forward and reverse barcodes. This many bases is removed from the
index 1 and index 2 of each read, concatenated, and used to demultiplex the
reads according to the provided map. (In our current pipeline, the indexes ARE
exactly this length.)

=item B<--troubleshoot_barcodes>

Try all four permutations of reverse-complementing and switching the 
concatenation order of the indexes. Whichever of the four permutations yields a
successful demux is used. Note: Performing these transformations on the indexes 
may coincidentally yield a barcode that seems to be correct, even though the 
overall demux is incorrect. 

=item B<-h>, B<--help>

Print help message and exit successfully.

=item B<--qsub-project> space, B<-qp> space

Indicate which qsub-project space should be used for all qsubmissions. The
default is jravel-lab.

=item B<--debug>, B<-d> {barcodes, demux, splitsamples, tagclean, dada2}

Runs one or more sections of the pipeline. To run multiple sections, type 
"--debug <section>" or "-d <section>" for each section. 

Section inputs:

Barcodes

1) raw reads and indexes
2) map

Demux

1) ./barcodes.fastq
2) map

Splitsamples

1) ./fwdSplit/seqs.fastq
2) ./revSplit/seqs.fastq

Tagclean

1) ./fwdSplit/split_by_sample_out/<sample_id>_*.fastq
2) ./revSplit/split_by_sample_out/<sample_id>_*.fastq

DADA2 (the sample_id is delimited by the first underscore character)

1) ./<sample_id>_*R1_tc.fastq
2) ./<sample_id>_*R2_tc.fastq

=item B<--verbose>

Prints every shell command to: <storage-dir>/<project>/<run>/qsub_stdout_logs>/illumina_dada2.pl.stdout

=item B<--dry-run>

Runs the pipeline without executing any of the shell commands. May be useful
combined with B<--verbose>. (Currently with B<--dry-run>, the pipeline may not 
progress far due to checkpoints that halt the pipeline if any step seems to 
fail.)

=item B<--email>

Notify by email when the job is finished. Does this by adding "-m ea" to the
outermost qsub call. Compatible with --qsub.

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

Most options will override the defaults shown above. The qsub options must be 
specified as a single string surrounded by double quotes ("), as shown below.

    part1.sh --qsub="-m ea -l excl=true" ...

=item B<--dada2>="options"

Overrides the default DADA2 parameters used at the MSL. The following options
are allowed:

 --dada2-truncLen-f, -for (defaults: V3V4: 225 | V4: 200 | ITS: 0)
 --dada2-truncLen-r, -rev (defaults: V3V4: 225 | V4: 200 | ITS: 0)
 --dada2-maxN (default: 0)
 --dada2-maxEE (defaults: V3V4: 2 | V4: 2 | ITS: 0)
 --dada2-truncQ (default: 2)
 --dada2-rmPhix (default: TRUE)
 --dada2-maxLen (default: Inf)
 --dada2-minLen (default: V3V4: 20 | V4: 20 | ITS: 50)
 --dada2-minQ (default: 0)

Please see https://rdrr.io/bioc/dada2/man/filterAndTrim.html for descriptions
of the parameters. The parameters should be given within double quotes as shown 
below:

part1.sh --dada2="--dada2-maxEE 5 --dada2-minQ 10" ...

=item B<--dada2-mem> memory

The amount of memory to request for the DADA2 qsub job. Use Sun Grid Engine
syntax: 10G, 500M, etc.

=back 

=cut

