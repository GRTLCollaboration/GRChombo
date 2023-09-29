#!/bin/bash

#! System to be used
#SBATCH -p skylake
#! Number of nodes to be used
#SBATCH -N 3
#! Number of MPI tasks (<=32*nodes)
#SBATCH --ntasks=90
#! Number of MPI tasks allocated to each node
#SBATCH --cpus-per-task=4
#! Wall clock time required
#SBATCH --time=12:00:00
#! Name of job
#SBATCH -J L8_long

srun ./Main_ScalarField3d_ch.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.ICPX2022.0.ex params.txt
