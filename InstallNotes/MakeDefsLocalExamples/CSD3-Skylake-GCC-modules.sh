#!/bin/bash
# To load the modules below call this file as an argument to `source`
module unload intel/bundles/complib/2017.4
# Load newer GCC version
module swap gcc gcc-7.2.0-gcc-4.8.5-pqn7o2k
module load openmpi-3.1.3-gcc-7.2.0-b5ihosk
module load hdf5-1.10.4-gcc-5.4.0-7zl2gou
module load intel/mkl/2019.3
