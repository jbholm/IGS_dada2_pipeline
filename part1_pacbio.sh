#!/bin/bash
set -e # exit when any command fails (returns non-0 exit status)
MY_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

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

# global params
source <(python3 $MY_DIR/lib/get_config.py modules pacbio_demux_denoise global)
# full-length is default
source <(python3 $MY_DIR/lib/get_config.py modules pacbio_demux_denoise full-length)

while [[ ! "$1" == "--" && "$#" != 0 ]]; do
  case "$1" in
    --qsub*)
        if [[ $1 =~ "--qsub=" ]]; then 
            QSUB_ARGS="${1#*=}"
            if [[ ! -n "$QSUB_ARGS" || ! $1 =~ "=" ]]; then  
                MSG="--qsub missing value."
                MSG+=" --qsub=\"\" and --qsub= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--qsub" ]]; then
            try_assign QSUB_ARGS "$1" "$2"
            shift 2
        fi
        ;;
    --args*)
        if [[ $1 =~ "--args=" ]]; then 
            OTHER_ARGS="${1#*=}"
            if [[ ! -n "$OTHER_ARGS" || ! $1 =~ "=" ]]; then  
                MSG="--args missing value."
                MSG+=" --args=\"\" and --args= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--args" ]]; then
            try_assign OTHER_ARGS "$1" "$2"
            shift 2
        fi
        ;;
    --pattern*) 
        if [[ $1 =~ "--pattern=" ]]; then 
            file_pattern="${1#*=}"
            if [[ ! -n "$file_pattern" || ! $1 =~ "=" ]]; then  
                MSG="--pattern missing value."
                MSG+=" --pattern=\"\" and --pattern= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--pattern" ]]; then
            try_assign file_pattern "$1" "$2"
            shift 2
        fi
        ;;
    --run|-r)
        if [[ $1 =~ "--run=" ]]; then 
            RUN="${1#*=}"
            if [[ ! -n "$RUN" || ! $1 =~ "=" ]]; then  
                MSG="--run missing value."
                MSG+=" --run=\"\" and --run= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--run" || $1 == "-r" ]]; then
            try_assign RUN "$1" "$2"
            shift 2
        fi
        ;;
    --run-storage*)
        if [[ $1 =~ "--run-storage=" ]]; then 
            SD="${1#*=}"
            if [[ ! -n "$SD" || ! $1 =~ "=" ]]; then  
                MSG="--run-storage missing value."
                MSG+=" --run-storage=\"\" and --run-storage= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--run-storage" ]]; then
            try_assign SD "$1" "$2"
            shift 2
        fi
        ;;
    --dada2_forward_primer*)
        if [[ $1 =~ "--dada2_forward_primer=" ]]; then 
            dada2_forward_primer="${1#*=}"
            if [[ ! -n "$dada2_forward_primer" || ! $1 =~ "=" ]]; then  
                MSG="--dada2_forward_primer missing value."
                MSG+=" --dada2_forward_primer=\"\" and --dada2_forward_primer= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--dada2_forward_primer" ]]; then
            try_assign dada2_forward_primer "$1" "$2"
            shift 2
        fi
        ;;
    --dada2_reverse_primer*)
        if [[ $1 =~ "--dada2_reverse_primer=" ]]; then 
            dada2_reverse_primer="${1#*=}"
            if [[ ! -n "$dada2_reverse_primer" || ! $1 =~ "=" ]]; then  
                MSG="--dada2_reverse_primer missing value."
                MSG+=" --dada2_reverse_primer=\"\" and --dada2_reverse_primer= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--dada2_reverse_primer" ]]; then
            try_assign dada2_reverse_primer "$1" "$2"
            shift 2
        fi
        ;;
    --dada2_min_length*)
        if [[ $1 =~ "--dada2_min_length=" ]]; then 
            dada2_min_length="${1#*=}"
            if [[ ! -n "$dada2_min_length" || ! $1 =~ "=" ]]; then  
                MSG="--dada2_min_length missing value."
                MSG+=" --dada2_min_length=\"\" and --dada2_min_length= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--dada2_min_length" ]]; then
            try_assign dada2_min_length "$1" "$2"
            shift 2
        fi
        ;;
    --dada2_max_length*)
        if [[ $1 =~ "--dada2_max_length=" ]]; then 
            dada2_max_length="${1#*=}"
            if [[ ! -n "$dada2_max_length" || ! $1 =~ "=" ]]; then  
                MSG="--dada2_max_length missing value."
                MSG+=" --dada2_max_length=\"\" and --dada2_max_length= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--dada2_max_length" ]]; then
            try_assign dada2_max_length "$1" "$2"
            dada2_max_length="--max_length=${dada2_max_length}"
            shift 2
        fi
        ;;
    --its|--ITS)
        source <(python3 $MY_DIR/lib/get_config.py modules pacbio_demux_denoise its)
        shift 1
        ;;
    --no-delete|--nodelete|--no_delete)
        NODELETE=$1
        shift 1
        ;;
    --email)
        EMAIL="-m ea"
        shift 1
        ;;
    --qsub-project*)
        if [[ $1 =~ "--qsub-project=" ]]; then 
            QP="${1#*=}"
            if [[ ! -n "$QP" || ! $1 =~ "=" ]]; then  
                MSG="--qsub-project missing value."
                MSG+=" --qsub-project=\"\" and --qsub-project= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--qsub-project" ]]; then
            try_assign QP "$1" "$2"
            shift 2
        fi
        ;;
    --help|-h)
        perldoc -F "${0}"
        exit 0
        ;;
    *) # preserve positional arguments even if they fall between other params
    # note that options and their operands will be separate elements of $PARAMS
      PARAMS="$PARAMS $1"
      shift 1
      ;;
  esac
done

# -r is mandatory
if [[ ! -n "$RUN" ]]; then
    stop "Run ID (-r) required."
fi

# Validate the storage directory
SD_DEFAULT=`cat "$MY_DIR/config.json" | \
    python3 -c "import sys, json; print(json.load(sys.stdin)['run_storage_path'])"`
SD_DEFAULT=$(readlink -f "$SD_DEFAULT")

if [[ ! -n "$SD" ]]; then
    SD=$SD_DEFAULT
elif [[ ! -d "$SD" ]]; then
    stop "$SD does not exist or is not a directory!\n" # A custom storage directory must exist
else
    SD="${SD%\/}" # Remove any trailing slash and get abs path
fi
SD=$(readlink -f "$SD/$RUN")
printf "%b" "WORKING DIRECTORY: $SD\n"

# get user to confirm overwrite if the run directory already exists
if [[ -d "$SD" ]]; then
    printf "$SD already exists!\n"
    select yn in 'Abort' 'Resume' 'Overwrite'; do 
        case "$yn" in 
            "Overwrite")  
                echo "Confirm complete overwrite? All existing progress will be lost."
                select yn2 in 'No' 'Yes'; do
                    case "$yn2" in
                        "Yes")
                        rm -rf $SD
                        mkdir -p $SD
                        break
                        ;;
                        "No")
                        stop
                        ;;
                    esac;
                done 
                break
                ;;
            "Abort")
                stop
                ;;
            "Resume")
                break
                ;;  
        esac; 
    done
else
    UMASK_OLD=`umask`
    if [ "${SD##$SD_DEFAULT}" != "$SD" ]; then  
        umask 0002 # if we're in the shared run directory, share the new directories
    fi
    mkdir -p "$SD"
fi
mkdir -p "$SD/qsub_error_logs/"
mkdir -p "$SD/qsub_stdout_logs/"
cd $SD

if [[ ! -n "$QP" ]]; then
    printf "qsub-project ID (--qsub-project) not provided. Using jravel-lab as default\n"
    QP="jravel-lab"
fi

if [[ -n "$NODELETE" ]]; then
    printf "%b\n" "Keeping intermediate files" 
fi


if [[ -n "$EMAIL" ]]; then
    QSUB_ARGS="$QSUB_ARGS $EMAIL"
fi

if [[ -n "$file_pattern" ]]; then
    file_pattern="--pattern=\"$file_pattern\""
fi

if [[ ! -z ${dada2_forward_primer+x} ]]; then
    if [[ ${dada2_forward_primer^^} =~ [^ACTGURYSWKMBVDHN] ]]; then
        stop "Unrecognized characters in --dada2_forward_primer"
    else
        dada2_forward_primer="--forward_primer=$dada2_forward_primer"
    fi
fi
if [[ ! -z ${dada2_reverse_primer+x} ]]; then
    if [[ ${dada2_reverse_primer^^} =~ [^ACTGURYSWKMBVDHN] ]]; then
        stop "Unrecognized characters in --dada2_reverse_primer"
    else
        dada2_reverse_primer="--reverse_primer=$dada2_reverse_primer"
    fi
fi
if [[ ! -z ${dada2_min_length+x} ]]; then
    if [[ ! ${dada2_min_length} =~ [0-9]+ ]]; then
        stop "--dada2_min_length must be integer"
    else
        dada2_min_length="--min_length=${dada2_min_length}"
    fi
fi
if [[ ! -z ${dada2_max_length+x} ]]; then
    if [[ ! ${dada2_max_length} =~ [0-9]+ ]]; then
        stop "--dada2_max_length must be integer"
    else
        dada2_max_length="--max_length=${dada2_max_length}"
    fi
fi

use () 
{ 
    eval `/usr/local/packages/usepackage/bin/usepackage -b $*` || true
}

#use sge > /dev/null 2>&1  || true
module load sge 2>/dev/null || true
export PYTHONPATH=""

#use r-3.6.0 > /dev/null 2>&1  || true
module load r/4.0.2 2>/dev/null || true

export LD_LIBRARY_PATH=/usr/lib64/:$LD_LIBRARY_PATH # for Rcpp libstdc++.so.6

# Begin log (will be continued by pacbio_dda2.R)
log="$SD/${RUN}_16S_pipeline_log.txt"

OPTSARR=("$PARAMS" "$NODELETE")
OPTS="${OPTSARR[*]}"
OPTS="$( echo "$OPTS" | awk '{$1=$1;print}' )" # remove excess spaces in above OPTSARR
R=`cat "$MY_DIR/config.json" | \
    python3 -sc "import sys, json; print(json.load(sys.stdin)['R'])"`
ARGS=("-l mem_free=128G" "-P" "$QP" "-V" "-N" "MSL_PACBIO" "-o ${SD}/qsub_stdout_logs/pacbio_dada2.R.stdout" "-e ${SD}/qsub_error_logs/pacbio_dada2.R.stderr" "$QSUB_ARGS" "${R}script" "$MY_DIR/pacbio_dada2.R" "$OPTS" "--wd" "$SD" "${file_pattern}" "${dada2_forward_primer}" "${dada2_reverse_primer}" "${dada2_min_length}" "${dada2_max_length}" "--rm.phix" "${OTHER_ARGS}")
CMD=()
for ARG in "${ARGS[@]}"; do # remove absent commands from above ARGS
    if [[ -n "$ARG" ]]; then
        CMD+=("$ARG")
    fi
done

QSUB_CMD=${CMD[@]}
QSUB_CMD="$QSUB_CMD"

EXECUTOR=`cat "$MY_DIR/config.json" | \
    python3 -sc "import sys, json; print(json.load(sys.stdin)['executor'])"`

printf "$ $EXECUTOR $QSUB_CMD\n"
printf "$ $EXECUTOR $QSUB_CMD\n" >> $log
$EXECUTOR $QSUB_CMD

: <<=cut
=pod

=head1 NAME

part1.sh

=head1 DESCRIPTION

The script can be launched from any location on the IGS server. It creates a
working directory in the global run location*. In the working directory, this
script:

1. Copies the demultiplexed raw reads into a ./demultiplexed/ directory. The
demultiplexed reads are found using the provided regex-like pattern.

2. Removes adapters from reads and writes new FASTQ's to ./tagcleaned/.

3. Filters reads on quality and writes new FASTQ's to ./filtered/.

4. Denoises reads using the dada2 pipeline with parameters specific to the
full-length 16S rRNA gene region.

A log file is written at ./<RUN>_16S_pipeline_log.txt. 

* See the pipeline's config.json for installation-specific variables.

=head1 SYNOPSIS

part1.sh <INPUT DIRECTORY> --pattern=<REGEX> -r <RUN> [<options>]

=head1 OPTIONS

=over

=item B<INPUT DIRECTORY>

Path to a directory where demultiplexed FASTQ's will be found.

=item B<--pattern> REGEX

The regex-like pattern that will be used to find demultiplexed FASTQ's and 
parse sample names. Use a named capture group in the form "(?<sample>...)" to 
denote where the sample name should be found. Give the pattern as a quoted 
string on the command line. See the config file for the default pattern.

=item B<-r> RUN

The run's name. The working directory will be named after the run, and the 
sample names will be appended to the run name, to uniquely identify samples
from different runs in Part 2.

=item b<--run-storage> PATH

Override the global run storage location. The working directory will be created
as a subdirectory of PATH.

=item B<--no-delete|--nodelete>

Don't delete intermediate files

=item B<--qsub> OPTIONS

A quoted string giving additional options to the qsub command. For example,
to start a run using less than the default RAM requirement, 
B<--qsub "-l mem_free=16G">

=item B<--email>

Appends "-m ea" to the qsub command so the current user is notified when the
job finishes.

=item B<--qsub-project> ID

Set the qsub project ID. Default: jravel-lab.

=head2 DADA2 PARAMETERS

=item B<--dada2_forward_primer>, <--dada2_reverse_primer>

Primer sequences to remove. Sequences without both primers are discarded. 
Reverse primer must be given as if read in the same direction as the forward 
primer.

=item B<--max_mismatch>=INTEGER

The number of mismatches to tolerate when matching reads to primer sequences.

=item B<--indels>

Allow insertions or deletions of bases when matching adapters

=item B<--dada2_max_length>=MAX_LENGTH

Remove reads with length longer than MAX_LENGTH. MAX_LENGTH is enforced before 
quality trimming and truncation.

=item B<--dada2_min_length>=MIN_LENGTH

Remove reads with length shorter than MIN_LENGTH. MIN_LENGTH is enforced after 
quality trimming and truncation.

=item B<--ITS|--its>

A synonym for B<--dada2_min_length=50>.

=back 

=cut
