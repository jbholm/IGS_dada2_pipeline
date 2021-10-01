#!/bin/sh
set -e

echo -e "Sample\tReads"
sed -n '/Median sequence length/,$p' - | sed '/Total number seqs written/q' - | head -n -2 | tail -n +2