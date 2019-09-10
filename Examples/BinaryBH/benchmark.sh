#!/bin/bash

outfile=$1

export OMP_NUM_THREADS=1
for binary in ../BH_0fast_native_vect ; do
  binary=$(basename $binary)
  for i in 1 2 3; do 
    mpiexec --report-bindings --bind-to core  --map-by core -np  64 ../${binary} ../params_expensive.txt
    printf "${binary} " >> $outfile
    grep 'Total simulation time' pout.0 >> $outfile
  done
done
