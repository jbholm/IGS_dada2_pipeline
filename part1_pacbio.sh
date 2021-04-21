#!/bin/bash
set -e # exit when any command fails (returns non-0 exit status)

use () 
{ 
    eval `/usr/local/packages/usepackage/bin/usepackage -b $*` || true
}
use r-3.6.0 > /dev/null 2>&1  || true
module load r/3.6.3 2>/dev/null || true

# use sge > /dev/null 2>&1  || true
# module load sge 2>/dev/null || true

MY_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
"$MY_DIR/pacbio_dada2.R" "$@"
# CMD="-cwd -w e -b y -P jravel-lab -q threaded.q -pe thread 4 -V -m ea -l mem_free=8G -N MSL_PACBIO ${MY_DIR}/pacbio_dada2.R"
# echo "$ qsub $CMD "
# echo "$@"
# echo "\n"
# qsub $CMD "$@"