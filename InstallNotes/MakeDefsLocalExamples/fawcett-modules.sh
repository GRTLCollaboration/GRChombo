#!/bin/bash
# To load the modules below call this file as an argument to `source`
module load gcc/7.5.0/gcc-7.5.0-4zdnknt # in default environment but load anyway in case of purge
module load intel-oneapi-compilers/2022.0.2/gcc-7.5.0-3clffac intel-oneapi-mkl/2022.0.2/gcc-7.5.0-rnpbvnf
module load hdf5/1.12.1/intel-oneapi-mpi-2021.5.1/gcc-7.5.0-34j4qx4
module load petsc/3.16.5/intel-oneapi-mpi-2021.5.1/gcc-7.5.0-67nf5rm # for AHFinder
