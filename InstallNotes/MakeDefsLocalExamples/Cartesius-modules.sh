module load fortran/intel/17.0.4
module load c/intel/17.0.4
module load gcc/5.2.0
module load hdf5/impi/intel/1.8.12
module load mpi/impi/17.0.2
module load mkl/17.0.4

export LD_LIBRARY_PATH=/hpc/sw/hdf5-1.8.12-intel-impi-par/lib:$LD_LIBRARY_PATH
source /opt/intel/parallel_studio_xe_2017_update4/mkl/bin/mklvars.sh intel64

