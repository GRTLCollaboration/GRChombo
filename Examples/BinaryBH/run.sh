#!?bin/bash

taskset -c -p 1 $$
./Main_BinaryBH3d.Linux.64.g++.gfortran.OPTHIGH.OPENMPCC.ex params_very_cheap.txt
