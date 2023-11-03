#!/bin/bash

#! System to be used
#SBATCH -p knl
#! Number of nodes to be used
#SBATCH -N 2
#! Number of MPI tasks (<=32*nodes)
#SBATCH --ntasks=8
#! Number of MPI tasks allocated to each node
#SBATCH --cpus-per-task=4
#! Wall clock time required
#SBATCH --time=02:00:00
#! Name of job
#SBATCH -J slurm-check-test

srun ./Main_ScalarField3d_ch.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.ICPX2022.0.ex params.txt
