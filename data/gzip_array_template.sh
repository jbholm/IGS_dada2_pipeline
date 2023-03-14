#!/bin/bash
#$ -l mem_free=1G

#$ -t 1-1756

infile="$(awk "NR==$SGE_TASK_ID" demuxed_seqs.lst)"
gzip -f9 $infile
