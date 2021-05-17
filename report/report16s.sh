#!/bin/bash
set -e # exit when any command fails (returns non-0 exit status)


export LD_LIBRARY_PATH=/usr/local/packages/gcc/lib64
export PATH=/usr/local/bin:$PATH

/home/jolim/share/miniconda3/envs/msl_reports/bin/python3.7 -s $(dirname "$0")/make_report.py "$@"

