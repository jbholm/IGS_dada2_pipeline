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


while [[ ! "$1" == "--" && "$#" != 0 ]]; do
  case "$1" in
    # --qsub*)
    #     QSUB_ARGS="${1#*=}"
    #     if [[ ! -n "$QSUB_ARGS" || ! $1 =~ "=" ]]; then  
    #         MSG="--qsub missing value."
    #         MSG+=" --qsub=\"\" and --qsub= are not accepted."
    #         stop "$MSG"
    #     fi
    #     shift 1
    #     ;;
    --pattern) 
        try_assign PATTERN "$1" "$2"
        shift 2
        ;;
    --run|-r)
        try_assign RUN "$1" "$2"
        shift 2
        ;;
    --run-storage)
        try_assign SD "$1" "$2"
        shift 2
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
if [[ ! -n "$SD" ]]; then
    SD=`cat "$MY_DIR/config.json" | \
    python3 -c "import sys, json; print(json.load(sys.stdin)['run_storage_path'])"`
elif [[ ! -d "$SD" ]]; then
    stop "$SD does not exist or is not a directory!\n" # A custom storage directory must exist
else
    SD="${SD%\/}" # Remove any trailing slash and get abs path
fi
SD=$(readlink -f "$SD/$RUN")
printf "%b" "WORKING DIRECTORY: $SD\n"
mkdir -p "$SD"
cd "$SD"

# get user to confirm overwrite if the run directory already exists
if [[ -d "$SD" ]]; then
    printf "$SD already exists! Overwrite?\n"
    select yn in 'Yes' 'No'; do 
        case "$yn" in 
            "Yes")  
                rm -rf "$SD/*"
                break
                ;; 
            "No")
                stop
                ;; 
        esac; 
    done
fi

use () 
{ 
    eval `/usr/local/packages/usepackage/bin/usepackage -b $*` || true
}
use r-3.6.0 > /dev/null 2>&1  || true
module load r/3.6.3 2>/dev/null || true

# use sge > /dev/null 2>&1  || true
# module load sge 2>/dev/null || true

# Begin log (will be continued by pacbio_dda2.R)
log="$SD/${RUN}_16S_pipeline_log.txt"

OPTSARR=("$PARAMS")
OPTS="${OPTSARR[*]}"
OPTS="$( echo "$OPTS" | awk '{$1=$1;print}' )"
ARGS=("$OPTS" "--run" "$RUN" "--pattern" "$PATTERN")
CMD=()
for ARG in "${ARGS[@]}"; do
    if [[ -n "$ARG" ]]; then
        CMD+=("$ARG")
    fi 
done

# printf "$ qsub ${CMD[*]}\n"
# printf "$ qsub ${CMD[*]}\n" >> $log
# qsub ${CMD[*]}


# Somehwat misleading because the --pattern is given without quotes, when a 
# bash command would require quotes. Yet, the pacbio_dada2.R does receive the
# operand properly
# printf "$ $MY_DIR/pacbio_dada2.R ${CMD[*]}\n"
# printf "$ $MY_DIR/pacbio_dada2.R ${CMD[*]}\n" >> $log
"$MY_DIR/pacbio_dada2.R" ${CMD[*]}
# CMD="-cwd -w e -b y -P jravel-lab -q threaded.q -pe thread 4 -V -m ea -l mem_free=8G -N MSL_PACBIO ${MY_DIR}/pacbio_dada2.R"
# echo "$ qsub $CMD "
# echo "$@"
# echo "\n"
# qsub $CMD "$@"