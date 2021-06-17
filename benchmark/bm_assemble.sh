#!/bin/bash

echo "Benchmarking the system matrix assembly."

program=assemble.out
coord_file=../inputs/grid2.coord.txt
topo_file=../inputs/grid2.topo.txt

echo "Coordinate file: ${coord_file}"
echo "Topology file:   ${topo_file}"

fp=[0-9]+\.[0-9]+[Ee][+-][0-9]+ # floating point in scientific notation
logfile=bm_${program}_$(date +%s).tsv

for i in {1..2}; do
    export OMP_NUM_THREADS=$i
    res=$(./${program} ${coord_file} ${topo_file} <<< $filename | grep "Elapsed time" | egrep -o $fp)
    printf "%i\t%s\n" $i "$res"
    printf "%i\t%s\n" $i "$res" >> $logfile
done

echo "Saved to $logfile"