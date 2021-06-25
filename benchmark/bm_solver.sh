#!/bin/bash

echo "Benchmarking the whole solver."

coord_file=../inputs/grid3.coord.txt
topo_file=../inputs/grid3.topo.txt

echo "Coordinate file: ${coord_file}"
echo "Topology file:   ${topo_file}"

fp=[0-9]+\.[0-9]+[Ee][+-][0-9]+ # floating point in scientific notation
logfile=bm_${program}_$(date +%s).tsv

for i in {1..8}; do
    export OMP_NUM_THREADS=$i
    res=$(../solver.out ${coord_file} ${topo_file} | grep "Elapsed time" | egrep -o $fp)
    printf "%i\t%s\n" $i "$res"
    printf "%i\t%s\n" $i "$res" >> $logfile
done

echo "Saved to $logfile"