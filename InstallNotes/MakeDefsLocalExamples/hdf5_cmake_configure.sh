#!/bin/zsh

# Script to configure hdf5 in parallel on Mac OS Catalina using gcc compilers
# see full instructions and obtain the files at 
# https://support.hdfgroup.org/HDF5/release/cmakebuild.html 
export MPI_DIR=/usr/local/Cellar/mpich/3.3.2
export CC=gcc-9
export CXX=g++-9
export FC=gfortran-9
export F90=gfortran-9
export F77=gfortran-9
export CFLAGS="-I${MPI_DIR}/include -L${MPI_DIR}/lib -lmpi"
export CXXFLAGS="-I${MPI_DIR}/include -L${MPI_DIR}/lib -lmpi"

./configure --enable-build-mode=production --enable-parallel --prefix=/usr/local/Cellar/hdf5-parallel 