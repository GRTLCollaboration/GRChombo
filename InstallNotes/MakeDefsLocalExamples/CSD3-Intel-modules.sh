#!/bin/bash
# To load the modules below call this file as an argument to `source`
module unload intel/bundles/complib/2017.4
# Load newer GCC version
module swap gcc gcc-7.2.0-gcc-4.8.5-pqn7o2k
# MR: I get MPI problems with 2019 so sticking with 2018 for now
module load intel/bundles/complib/2018.4
module load hdf5-1.10.4-intel-17.0.4-swn7n43
