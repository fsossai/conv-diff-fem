#!/bin/bash

if (($# < 1)); then
    echo "Usage: $0 <program>"
    exit 1
fi

if [ ! -f $1 ]; then
    echo "ERROR: '$1' file not found"
    exit 1
fi

program=$1
filename="../inputs/mat1500.txt"
fp=[0-9]+\.[0-9]+[Ee][+-][0-9]+ # floating point in scientific notation
logfile=bm_${program}_$(date +%s).tsv

for i in {1..4}; do
    export OMP_NUM_THREADS=$i
    res=$(./${program} <<< $filename | grep "Elapsed time" | egrep -o $fp)
    printf "%i\t%s\n" $i "$res"
    printf "%i\t%s\n" $i "$res" >> $logfile
done

echo Saved to $logfile