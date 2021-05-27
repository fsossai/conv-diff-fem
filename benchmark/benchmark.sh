#!/bin/bash

program=amxpby
iterations=100
fp=[0-9]+\.[0-9]+[Ee][+-][0-9]+ # floating point in scientific notation
logfile=bm_${program}_$(date +%s).tsv

for i in {1..8}; do
    export OMP_NUM_THREADS=$i
    res=$(./${program}.out <<< $iterations | grep GFlops | egrep -o $fp)
    printf "%i\t%s\n" $i "$res"
    printf "%i\t%s\n" $i "$res" > $logfile
done

echo Saved to $logfile