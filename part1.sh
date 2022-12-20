#!/bin/bash
set -e # exit when any command fails (returns non-0 exit status)
# ensure the current working directory exists. (`module load sge` will fail if
# working directory is stale file handle.)
cd ~ && cd - 1>/dev/null

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
    --help|-h)
        perldoc -F "${0}"
        exit 0
        ;;
    -r|--run-ID)
        try_assign RUN "$1" "$2"
        shift 2
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
        else
            PARAMS="$PARAMS $1"
            shift 1
        fi
        ;;
    --1Step)
        ONESTEP=$1
        shift 1
        ;;
    --verbose)
        VERBOSE=$1
        shift 1
        ;;
    --dry-run)
        DRY_RUN=$1
        shift 1
        ;;
    --debug|-d)
        if [[ "$2" =~ ^- || ! -n "$2" ]]; then
            stop "--debug missing its value. Unable to continue."
        else 
            DBG="$DBG --debug $2"
            shift 2
        fi
        ;;
    --noskip|--no-skip)
        NOSKIP=$1
        shift 1
        ;;
    --no-delete|--nodelete)
        NODELETE=$1
        shift 1
        ;;
    --email)
        EMAIL="-m ea"
        shift 1
        ;;
    -qp|--qsub-project*)
        if [[ $1 =~ "--qsub-project=" ]]; then 
            QP="${1#*=}"
            if [[ ! -n "$QP" || ! $1 =~ "=" ]]; then  
                MSG="--qsub-project missing value."
                MSG+=" --qsub-project=\"\" and --qsub-project= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--qsub-project" || $1 == "-qp" ]]; then
            try_assign QP "$1" "$2"
            shift 2
        else
            PARAMS="$PARAMS $1"
            shift 1
        fi
        ;;
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
        else
            PARAMS="$PARAMS $1"
            shift 1
        fi
        ;;
    -i)
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
    -m|--map)
        try_assign MAP "$1" "$2"
        shift 2
        ;;
    --bclen*)
        if [[ $1 =~ "--bclen=" ]]; then 
            BCLENGTH="${1#*=}"
            if [[ ! -n "$BCLENGTH" || ! $1 =~ "=" ]]; then  
                MSG="--bclen missing value."
                MSG+=" --bclen=\"\" and --bclen= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--bclen" ]]; then
            try_assign BCLENGTH "$1" "$2"
            shift 2
        else
            PARAMS="$PARAMS $1"
            shift 1
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
    --fwd_primer*)
        if [[ $1 =~ "--fwd_primer=" ]]; then 
            FWD_PRIMER="${1#*=}"
            shift 1
        elif [[ $1 == "--fwd_primer" ]]; then
            try_assign FWD_PRIMER "$1" "$2"
            shift 2
        else
            PARAMS="$PARAMS $1"
            shift 1
        fi
        ;;
    --rev_primer*)
        if [[ $1 =~ "--rev_primer=" ]]; then 
            REV_PRIMER="${1#*=}"
            shift 1
        elif [[ $1 == "--rev_primer" ]]; then
            try_assign REV_PRIMER "$1" "$2"
            shift 2
        else
            PARAMS="$PARAMS $1"
            shift 1
        fi
        ;;
    --trim-maxlength*)
        if [[ $1 =~ "--trim-maxlength=" ]]; then 
            TRIM_MAXLENGTH="${1#*=}"
            if [[ ! -n "$TRIM_MAXLENGTH" || ! $1 =~ "=" ]]; then  
                MSG="--bclen missing value."
                MSG+=" --bclen=\"\" and --bclen= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--trim-maxlength" ]]; then
            try_assign TRIM_MAXLENGTH "$1" "$2"
            shift 2
        else
            PARAMS="$PARAMS $1"
            shift 1
        fi
        ;;
    --dada2-truncQ*)
        if [[ $1 =~ "--dada2-truncQ=" ]]; then 
            DADA2_TRUNCQ="${1#*=}"
            if [[ ! -n "$DADA2_TRUNCQ" || ! $1 =~ "=" ]]; then  
                MSG="--dada2-truncQ missing value."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--dada2-truncQ" ]]; then
            try_assign DADA2_TRUNCQ "$1" "$2"
            shift 2
        else
            PARAMS="$PARAMS $1"
            shift 1
        fi
        ;;
    --dada2-truncLen-f*)
        if [[ $1 =~ "--dada2-truncLen-f=" ]]; then 
            TRIM_LENGTH_FWD="${1#*=}"
            if [[ ! -n "$TRIM_LENGTH_FWD" || ! $1 =~ "=" ]]; then  
                MSG="--dada2-truncLen-f missing value."
                MSG+=" --dada2-truncLen-f=\"\" and --dada2-truncLen-f= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--dada2-truncLen-f" ]]; then
            try_assign TRIM_LENGTH_FWD "$1" "$2"
            shift 2
        else
            PARAMS="$PARAMS $1"
            shift 1
        fi
        ;;
    --dada2-truncLen-r*)
        if [[ $1 =~ "--dada2-truncLen-r=" ]]; then 
            TRIM_LENGTH_REV="${1#*=}"
            if [[ ! -n "$TRIM_LENGTH_REV" || ! $1 =~ "=" ]]; then  
                MSG="--dada2-truncLen-r missing value."
                MSG+=" --dada2-truncLen-r=\"\" and --dada2-truncLen-r= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--dada2-truncLen-r" ]]; then
            try_assign TRIM_LENGTH_REV "$1" "$2"
            shift 2
        else
            PARAMS="$PARAMS $1"
            shift 1
        fi
        ;;
    --amplicon_length*)
        if [[ $1 =~ "--amplicon_length=" ]]; then 
            AMPLICON_LENGTH="${1#*=}"
            if [[ ! -n "$AMPLICON_LENGTH" || ! $1 =~ "=" ]]; then  
                MSG="--amplicon_length missing value."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--amplicon_length" ]]; then
            try_assign AMPLICON_LENGTH "$1" "$2"
            shift 2
        else
            PARAMS="$PARAMS $1"
            shift 1
        fi
        ;;
    --dada2-minLen*)
        if [[ $1 =~ "--dada2-minLen=" ]]; then 
            DADA2_MINLEN="${1#*=}"
            if [[ ! -n "$DADA2_MINLEN" || ! $1 =~ "=" ]]; then  
                MSG="--dada2-minLen missing value."
                MSG+=" --dada2-minLen=\"\" and --dada2-minLen= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--dada2-minLen" ]]; then
            try_assign DADA2_MINLEN "$1" "$2"
            shift 2
        else
            PARAMS="$PARAMS $1"
            shift 1
        fi
        ;;
    --dada2-minQ*)
        if [[ $1 =~ "--dada2-minQ=" ]]; then 
            DADA2_MINQ="${1#*=}"
            if [[ ! -n "$DADA2_MINQ" || ! $1 =~ "=" ]]; then  
                MSG="--dada2-minQ missing value."
                MSG+=" --dada2-minQ=\"\" and --dada2-minQ= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--dada2-minQ" ]]; then
            try_assign DADA2_MINQ "$1" "$2"
            shift 2
        else
            PARAMS="$PARAMS $1"
            shift 1
        fi
        ;;
    --dada2-maxEE*)
        if [[ $1 =~ "--dada2-maxEE=" ]]; then 
            DADA2_MAXEE="${1#*=}"
            if [[ ! -n "$DADA2_MAXEE" || ! $1 =~ "=" ]]; then  
                MSG="--dada2-maxEE missing value."
                MSG+=" --dada2-maxEE=\"\" and --dada2-maxEE= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--dada2-maxEE" ]]; then
            try_assign DADA2_MAXEE "$1" "$2"
            shift 2
        else
            PARAMS="$PARAMS $1"
            shift 1
        fi
        ;;
    --dada2-rmPhix)
        DADA2_RMPHIX=$1
        shift 1
        ;;
    --no-dada2-rmPhix)
        DADA2_RMPHIX=$1
        shift 1
        ;;
    --dada2-mem*)
        if [[ $1 =~ "--dada2-mem=" ]]; then 
            DADA2MEM="${1#*=}"
            if [[ ! -n "$DADA2MEM" || ! $1 =~ "=" ]]; then  
                MSG="--dada2-mem missing value."
                MSG+=" --dada2-mem=\"\" and --dada2-mem= are not accepted."
                stop "$MSG"
            fi
            shift 1
        elif [[ $1 == "--dada2-mem" ]]; then
            try_assign DADA2MEM "$1" "$2"
            shift 2
        else
            PARAMS="$PARAMS $1"
            shift 1
        fi
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
if [[ -z ${MAP+x} ]]; then
    stop "\tPlease provide a full path to the run mapping file (-m)"
elif [[ ! -e "$MAP" ]]; then
    stop "\tMapping file (-m) does not exist or is inaccessible."
fi
printf "MAPPING FILE: $MAP\n"

if [[ ! -z ${BCLENGTH+x} ]]; then
    if [[ "$BCLENGTH" =~ [0-9]+ ]]; then
        BCLENGTH="--bclen=$BCLENGTH"
    else
        stop "--bclen must be a positive integer"
    fi
fi

# -r is mandatory
if [[ -z ${RUN+x} ]]; then
    stop "Run ID (-r) required."
fi
# if [[ $RUN =~ "[^a-zA-Z0-9\.]" ]]; then
#     stop "Run can only contain alphanumeric and '.' characters."
# fi

# Validate the variable region
[[ $ONESTEP ]] && onestep_key="onestep" || onestep_key="twostep"
supported_regions=( `cat "$MY_DIR/config.json" | \
    python3 -sc "import sys, json; \
    regions = json.load(sys.stdin)['part1 params']['$onestep_key'].keys(); \
    print(' '.join(regions))"` )
if [[ -z ${VAR+x} ]]; then
    VAR="V3V4"
else
    for reg in "${supported_regions[@]}"; do
        if [[ ${VAR^^} == $reg ]]; then
            supported="t"
        fi
    done
    if [[ ! -n $supported ]]; then
        if [[ ! -n $FWD_PRIMER || ! -n $REV_PRIMER || ! ( -n $TRIM_LENGTH_FWD && -n $TRIM_LENGTH_REV || -n $AMPLICON_LENGTH  ) ]]; then
            MSG="Variable region was '$VAR' but only the following are supported (case-insensitive):\n"
            MSG+="${supported_regions[*]}\n"
            MSG+="For unsupported regions, the pipeline requires:\n--fwd_primer\n--rev_primer\n"
            MSG+="either --fwd_trim_length and --rev_trim_length, or --amplicon_length.\n"
            stop "$MSG"
        fi
    fi
fi
printf "VARIABLE REGION: $VAR\n"
if [[ ! -z ${FWD_PRIMER+x} ]]; then
    if [[ $FWD_PRIMER =~ [^ACTGactg] ]]; then
        stop "Unrecognized characters in --fwd_primer"
    else
        FWD_PRIMER="--fwd_primer=$FWD_PRIMER"
    fi
fi
if [[ ! -z ${REV_PRIMER+x} ]]; then
    if [[ $REV_PRIMER =~ [^ACTGactg] ]]; then
        stop "Unrecognized characters in --rev_primer"
    else
        REV_PRIMER="--rev_primer=$REV_PRIMER"
    fi
fi
if [[ ! -z ${TRIM_LENGTH_FWD+x} ]]; then
    if [[ "$TRIM_LENGTH_FWD" =~ [0-9]+ ]]; then
        TRIM_LENGTH_FWD="--dada2-truncLen-f=$TRIM_LENGTH_FWD"
    else
        stop "--dada2-truncLen-f must be a positive integer"
    fi
fi
if [[ ! -z ${TRIM_LENGTH_REV+x} ]]; then
    if [[ "$TRIM_LENGTH_REV" =~ [0-9]+ ]]; then
        TRIM_LENGTH_REV="--dada2-truncLen-r=$TRIM_LENGTH_REV"
    else
        stop "--dada2-truncLen-r must be a positive integer"
    fi
fi
if [[ ! -z ${AMPLICON_LENGTH+x} ]]; then
    if [[ "$AMPLICON_LENGTH" =~ [0-9]+ ]]; then
        AMPLICON_LENGTH="--amplicon_length=$AMPLICON_LENGTH"
    else
        stop "--amplicon_length must be a positive integer"
    fi
fi
if [[ ! -z ${TRIM_MAXLENGTH+x} ]]; then
    if [[ "$TRIM_MAXLENGTH" =~ [0-9]+ ]]; then
        TRIM_MAXLENGTH="--trim-maxlength=$TRIM_MAXLENGTH"
    else
        stop "--trim-maxlength must be a positive integer"
    fi
fi

if [[ ! -z ${DADA2_TRUNCQ+x} ]]; then
    if [[ "$DADA2_TRUNCQ" =~ [0-9]+ ]]; then
        DADA2_TRUNCQ="--dada2-truncQ=$DADA2_TRUNCQ"
    else
        stop "--dada2-truncQ must be an integer"
    fi
fi

if [[ ! -z ${DADA2_MINLEN+x} ]]; then
    if [[ "$DADA2_MINLEN" =~ [0-9]+ ]]; then
        DADA2_MINLEN="--dada2-minLen=$DADA2_MINLEN"
    else
        stop "--dada2-minLen must be an integer"
    fi
fi
if [[ ! -z ${DADA2_MINQ+x} ]]; then
    if [[ "$DADA2_MINQ" =~ [0-9]+ ]]; then
        DADA2_MINQ="--dada2-minQ=$DADA2_MINQ"
    else
        stop "--dada2-minQ must be an integer"
    fi
fi
if [[ ! -z ${DADA2_MAXEE+x} ]]; then
    if [[ "$DADA2_MAXEE" =~ [0-9\.]+ ]]; then
        DADA2_MAXEE="--dada2-maxEE=$DADA2_MAXEE"
    else
        stop "--dada2-maxEE must be a number"
    fi
fi

SD_DEFAULT=`cat "$MY_DIR/config.json" | \
    python3 -c "import sys, json; print(json.load(sys.stdin)['run_storage_path'])"`
SD_DEFAULT=$(readlink -f "$SD_DEFAULT")
# Validate the storage directory
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
fi

if [[ ! -z ${DADA2MEM+x} ]]; then
    DADA2MEM="--dada2-mem=${DADA2MEM}"
fi

UMASK_OLD=`umask`
if [ "${SD##$SD_DEFAULT}" != "$SD" ]; then  
    umask 0002 # if we're in the shared run directory, share the new directories
fi
mkdir -p "$SD/qsub_error_logs/"
mkdir -p "$SD/qsub_stdout_logs/"

if [[ ! -n "$QP" ]]; then
    printf "qsub-project ID (--qp) not provided. Using jravel-lab as default\n"
    QP=""
else
    QP="-P $QP"
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

if [[ -n "$NODELETE" ]]; then
    printf "%b\n" "Keeping intermediate files" 
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
export PYTHONPATH=""
# if [[ -e "/usr/local/packages/miniconda3/etc/profile.d/conda.sh" ]]; then
# all this is safe because the environment changes won't persist outside of this
# shell!
# conda_init=`cat "$MY_DIR/config.json" | \
#     python3 -c "import sys, json; print(json.load(sys.stdin)['conda_init'])"`
# source "$conda_init" 2>/dev/null
# qiime_env=`cat "$MY_DIR/config.json" | \
#     python3 -c "import sys, json; print(json.load(sys.stdin)['qiime_env'])"`
# conda activate "$qiime_env"
# export PATH=/usr/local/packages/python-2.7/bin:$PATH
# export LD_LIBRARY_PATH=/usr/local/packages/python-2.7/lib:/usr/local/packages/gcc/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/lib64/:$LD_LIBRARY_PATH # for Rcpp libstdc++.so.6

module load r/4.0.2 2>/dev/null || true

# # Begin log (will be continued by illumina_dada2.pl)
log="$SD/${RUN}_16S_pipeline_log.txt"

# Remove extra spaces caused by joining empty arguments with a whitespace
OPTSARR=("$PARAMS" "$DBG" "$NOSKIP" "$NODELETE" "$VERBOSE" "$BCLENGTH" "$TROUBLESHOOT_BARCODES" "$ONESTEP" "$FWD_PRIMER" "$REV_PRIMER" "$TRIM_MAXLENGTH" "$DADA2_TRUNCQ" "$TRIM_LENGTH_FWD" "$TRIM_LENGTH_REV" "$AMPLICON_LENGTH" "$DADA2_MINLEN" "$DADA2_MINQ" "$DADA2_MAXEE" "$DADA2_RMPHIX" "$DADA2MEM" "$DRY_RUN")
OPTS="${OPTSARR[*]}"
OPTS="$( echo "$OPTS" | awk '{$1=$1;print}' )"

ARGS=("-l mem_free=8G" "-V" "$QP" "-N" "MSL_$RUN" "-o ${SD}/qsub_stdout_logs/illumina_dada2.pl.stdout" "-e ${SD}/qsub_error_logs/illumina_dada2.pl.stderr" "$QSUB_ARGS" "${MY_DIR}/illumina_dada2.pl" "$INPUT" "-wd" "$SD" "-v" "$VAR" "-m" "$MAP" "$OPTS")
CMD=()
for ARG in "${ARGS[@]}"; do
    if [[ -n "$ARG" ]]; then
        CMD+=("$ARG")
    fi 
done

EXECUTOR=`cat "$MY_DIR/config.json" | \
    python3 -sc "import sys, json; print(json.load(sys.stdin)['executor'])"`

printf "$ $EXECUTOR ${CMD[*]}\n"
printf "$ $EXECUTOR ${CMD[*]}\n" >> $log
$EXECUTOR ${CMD[*]}

: <<=cut
=pod

=head1 NAME

part1.sh


=head1 SYNOPSIS

part1.sh (-i <input directory> | -r1 <fwd reads> -r2 <rev reads> [-i1 <index 1> -i2 <index 2>]) -r <run> -m <map> [<options>]

=head1 DESCRIPTION

This is a wrapper to the MSL 16S pipeline for ILLUMINA runs. It can be launched 
from any location on the IGS server. When running the pipeline in full, the
essential inputs are a set of raw Illumina sequencing files (usually R1, R2, I1, 
and I2), a QIIME-formatted mapping file, and a run ID that will be used to name
the run directory.

The steps of the pipeline are as follows:

=over 

=item 1. Extracts barcodes from the raw files. Required inputs:

=over

=item Z<>* raw reads and indexes

=item Z<>* map

=back

=item 2. Demultiplexes the raw reads into fastq files containing reads specific to 
this project. Required inputs:

=over

=item Z<>* ./barcodes.fastq

=item Z<>* map

=back

=item 3. Produces individual .fastq files for each sample listed in the mapping file.
Required inputs:

=over

=item Z<>* ./fwdSplit/seqs.fastq

=item Z<>* ./revSplit/seqs.fastq

=back

=item 4. Performs tag-cleaning of each sample-specific file. Required inputs:

=over

=item Z<>* ./fwdSplit/split_by_sample_out/<sample_id>_*R1.fastq

=item Z<>* ./revSplit/split_by_sample_out/<sample_id>_*R2.fastq

=back

=item 5. Runs the forward and reverse reads through the dada2 pipeline for the 
V3V4 16S rRNA gene region. Alternatively, analysis of the V4 or ITS region may
be specified.

=over

=item Z<>* ./<sample_id>_*R1_tc.fastq

=item Z<>* ./<sample_id>_*R2_tc.fastq

=back

=back

Most pipeline products are stored in a directory named after the run. By default,
run directories are stored in /local/projects-t3/MSL/runs/. A log file is 
written at ./<RUN>_16S_pipeline_log.txt. Barcodes, raw libraries, and QC'd
reads are deleted if the DADA2 count table and stats exist when the pipeline
terminates.

=head1 OPTIONS

=head2 GENERAL

=over

=item B<-h>, B<--help>

Print help message and exit.

=item B<--run-ID>, B<-r> name

Create the run folder with this name.

=item B<--run-storage> path

Indicate an existing directory in which to place the run directory. The default
path is in the pipeline configuration file.

=item B<--1Step>

Use this flag if the data are prepared by 1-Step PCR (only r1 & r2 raw files
available)

=item B<--verbose>

Prints every shell command to the log file and to <run_directory>/qsub_stdout_logs/illumina_dada2.pl.stdout

=item B<--dry-run>

Runs the pipeline without executing any of the shell commands. May be useful
combined with B<--verbose>. (Currently with B<--dry-run>, the pipeline may not 
progress far due to checkpoints that halt the pipeline if any step seems to 
fail.)

=item B<--debug>, B<-d> {barcodes, demux, splitsamples, tagclean, dada2}

Runs the specified section of the pipeline. Multiple --debug options can be given
to run multiple consecutive parts of the pipeline, provided that the input to
the earliest requested step is present. Any non-consecutive steps will be 
ignored.
default is jravel-lab.

=item B<--no-skip|--noskip>

Will not check for any possible existing output files when deciding whether to 
run a section of the pipeline.

=item B<--no-delete|--nodelete>

Don't delete intermediate files

=item B<--email>

Notify by email when the job is finished. Does this by adding "-m ea" to the
outermost qsub call. Compatible with --qsub.

=item B<--qsub-project>, B<-qp> space

Indicate which qsub-project space should be used for all qsubmissions. The

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

=back

=head2 INPUT

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

=item B<--map>, B<-m> file

The full path to the Qiime-formatted mapping file.

=back

=head2 BARCODE EXTRACTION AND DEMULTIPLEXING

=over

=item B<--bclen> LENGTH

Manually specify the length of forward and reverse barcodes. This many bases is 
removed from the index 1 and index 2 of each read, concatenated, and used to 
demultiplex the reads according to the provided map.

In our current sequencing configuration, the indexes ARE exactly this length.
By default, the pipeline sets the barcode length equal to the length of the
first index.

=item B<--troubleshoot_barcodes>

Try all four permutations of reverse-complementing and switching the 
concatenation order of the indexes. Whichever of the four permutations yields a
successful demux is used. Note: Performing these transformations on the indexes 
may coincidentally yield a barcode that seems to be correct, even though the 
overall demux is incorrect. 

=back

=head2 TRIMMING, FILTERING, AND DENOISING

=over

This pipeline stores configurations for various hypervariable regions in its 
config file, located at <pipeline_path>/config.json. Some configurations are 
specific to one-step PCR library prep too. When the targeted variable region is 
not described in the config file, use the following options to specify the 
processing parameters. These options will also override any parameters specified 
by the config file.

=item B<--var-reg>, B<-v> {V3V4, V4, ITS}

The targeted variable region. V3V4 is default.

=item B<--fwd_primer>="SEQUENCE"

The primer(s) to trim from the forward reads. This value is passed on to 
bbduk.sh's "literal" option, so IUPAC ambiguity codes and comma-separated 
sequences are allowed.

=item B<--rev_primer>="SEQUENCE"

The primer(s) to trim from the reverse reads. See the note about B<--fwd_primer>.

=item B<--trim-maxlength>=LENGTH

After trimming primers, filter out reads longer than LENGTH. This value is
passed on to bbduk.sh's maxlen option. When primer trimming could yield various 
output read lengths, set this parameter to the maximum expected output read 
length.

=item B<--dada2-truncQ>=SCORE

Trim reads at the first instance of a quality score less than or equal to SCORE.

=item B<--dada2-truncLen-f>=LENGTH

Forward reads will be trimmed to LENGTH. Reads shorter than LENGTH will be
removed.

=item B<--dada2-truncLen-r>=LENGTH

Reverse reads will be trimmed to LENGTH. Reads shorter than LENGTH will be 
removed.

=item B<--amplicon_length>=LENGTH

The expected amplicon length. When determining the ideal trimming lengths, any
combination that prevents mate pair merging (due to excessive trimming) will be 
avoided. Amplicon length is not needed when B<--dada2-truncLen-f> and 
B<--dada2-truncLen-R> are both given.

=item B<--dada2-minLen>=LENGTH

Remove reads with length less than LENGTH. This filter is enforced after 
trimming.

=item B<--dada2-minQ>=SCORE

After trimming, reads contain a quality score less than SCORE will be discarded.

=item B<--dada2-maxEE>=VALUE

After trimming, reads with higher than VALUE "expected errors" will be 
discarded. 

=item B<--dada2-rmPhix>, B<--no-dada2-rmPhix>

Remove reads aligning to the Phi.X genome. (Give --no-dada2-rmPhix to override
a default value in the config file.)

=item B<--dada2-mem> memory

The amount of memory to request for the DADA2 qsub job. Use Sun Grid Engine
syntax: 10G, 500M, etc.

=back 

=cut

