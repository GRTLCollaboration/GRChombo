#!/bin/bash

#! System to be used
#SBATCH -p knl
#! Number of nodes to be used
#SBATCH -N 3
#! Number of MPI tasks (<=32*nodes)
#SBATCH --ntasks=90
>>>>>>> 12b652c (Pulling from Fawcett to resolve first slice load error)
#! Number of MPI tasks allocated to each node
#! Multiplication of N used above
#SBATCH --ntasks-per-node=8
>>>>>>> 9eada3d (Reset to previous commit because of issue reading in IC files, and fixed memory issue by de-allocating h vector?)
#! Wall clock time required
#SBATCH --time=02:00:00
#! Name of job
#SBATCH -J check000
>>>>>>> 9eada3d (Reset to previous commit because of issue reading in IC files, and fixed memory issue by de-allocating h vector?)

srun ./Main_ScalarField3d_ch.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.ICPX2022.0.ex params.txt
