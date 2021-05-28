#!/bin/bash

program=bicgstab
filename="../inputs/mat1500.txt"
fp=[0-9]+\.[0-9]+[Ee][+-][0-9]+ # floating point in scientific notation
logfile=bm_${program}_$(date +%s).tsv

for i in {1..4}; do
    export OMP_NUM_THREADS=$i
    res=$(./${program}.out <<< $filename | grep "Elapsed time" | egrep -o $fp)
    printf "%i\t%s\n" $i "$res"
    printf "%i\t%s\n" $i "$res" >> $logfile
done

echo Saved to $logfile