#!/bin/bash

filename=../inputs/Mgrid2.txt

if [ ! -f $filename ]; then
    echo "ERROR: '$filename' file not found"
    exit 1
fi

program=bicgstab.out
max_it=500

fp=[0-9]+\.[0-9]+[Ee][+-][0-9]+ # floating point in scientific notation
logfile=bm_${program}_$(date +%s).tsv

for i in {1..2}; do
    export OMP_NUM_THREADS=$i
    res=$(./${program} $filename $max_it | grep "Elapsed time" | egrep -o $fp)
    printf "%i\t%s\n" $i "$res"
    printf "%i\t%s\n" $i "$res" >> $logfile
done

echo "Saved to $logfile."