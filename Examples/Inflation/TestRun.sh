#!/bin/bash -l

#! How many nodes
#SBATCH --nodes 2
### NB cosma8 has 128 cores per node so product of these = 128
#SBATCH --ntasks-per-node=64 ## Total tasks will be this * number of nodes
#SBATCH --cpus-per-task=2 ## Equal to number of OMP threads per task (set below)
#! How many CPUs per tasks should be allocated (for OpenMP threading)
#SBATCH -J InflationEx
#SBATCH -o out_%J.out
#SBATCH -e err_%J.err
#SBATCH -p cosma8
#SBATCH -A dp092
#SBATCH --exclusive
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE                          # notifications for job done & fail
#SBATCH --mail-user=a.waeming@qmul.ac.uk

module purge
#### Cosma 7
##module load intel_comp/2019
##module load intel_mpi/2019
##module load parallel_hdf5/1.10.3
#### Cosma 8
module load gnu_comp/7.3.0
module load intel_comp/2022.1.2
module load compiler/2022.0.2 mkl/2022.0.2 mpi/2021.5.1
module load parallel_hdf5/1.12.0

#!Print info
module list
pwd
date

#! Are you using OpenMP?
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#! Number of nodes and tasks per node allocated by SLURM (do not change):
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$[${numnodes}*${mpi_tasks_per_node}]

#! Full path to application executable: 
##application="/cosma/home/dp092/dc-waem3/ModGrav_Reheating/GRFolres/Examples/Reheating/Main_KerrBH4dST3d.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.ex"
application="/cosma/home/dp092/dc-waem3/GRChomboForked/GRChombo/Examples/Inflation/Main_Inflation3d_ch.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.COSMA8.Intel2022.ex"
#! Run options for the application:
options="params.txt"

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"

# Run the program

mpirun -ppn $mpi_tasks_per_node -np $SLURM_NTASKS $application $options
