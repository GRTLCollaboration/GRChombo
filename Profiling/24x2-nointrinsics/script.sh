#!/bin/bash

export OMP_NUM_THREADS=2
export OMP_PROC_BIND=TRUE
export OMP_PLACES=cores
mpiexec -n 24 -cpus-per-proc 2 ../Main_BinaryBH3d_ch.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.OPENMPCC.GCC.11.0.malloc.nointrinsics.ex ../params.txt
