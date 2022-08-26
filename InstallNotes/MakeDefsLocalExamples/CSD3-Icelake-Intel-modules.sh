#!/bin/bash
# To load the modules below call this file as an argument to `source`

# Check that we are running on an Icelake (either login or compute) node
if [[ $(hostname) != *"-q-"* ]]; then
	echo "These modules will only work on an Icelake (Rocky Linux 8) node." 1>&2
	echo "See https://docs.hpc.cam.ac.uk/hpc/user-guide/connecting.html"
	echo "Not loading anything..." 1>&2
	return 1
fi

# This should already be loaded but load anyway in case modules have been purged
module load rhel8/default-icl # loads Intel compiler and MPI

module load intel-oneapi-mkl/2022.1.0/intel/mngj3ad6
module load hdf5/1.10.8/intel/intel-oneapi-mpi/h75adcal

# Uncomment if using AHFinder
#module load petsc/3.17-icl

# Uncomment if using TwoPunctures
#module load gsl/2.7.1/gcc/cjb5f4ui
