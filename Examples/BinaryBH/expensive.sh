#!/bin/bash
#SBATCH --account m1411
#SBATCH --qos=debug
#SBATCH --time=30
#SBATCH --nodes=8
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=4
#SBATCH --constraint=haswell

export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NUM_THREADS=4

srun --cpu-bind=cores Main_BinaryBH3d.Linux.64.CC.ftn.DEBUG.OPT.MPI.ex params_expensive.txt
