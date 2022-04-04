#!/bin/bash

#Load modules
module --force purge
#export MPICH_CXX=icpc
module load modules/1.59-20220201
#module load intel-oneapi-compilers/2022.0.1 intel-oneapi-mkl/2022.0.1 intel-oneapi-mpi/2021.5.0 fftw/3.3.10-mpi hdf5/1.12.1-mpi
module load intel-oneapi-compilers intel-oneapi-mkl openmpi/4.0.7 fftw/3.3.10-mpi hdf5/1.8.22-mpi

# Run configuration script and save to logfile
#./configure.py --prob=02_test_turb_cooling -fft -mpi --cxx icpc -hdf5 -h5double > configure.log 
#./configure.py --prob=04_test_turb_sn_cooling -fft -mpi -hdf5 -h5double > configure.log
./configure.py --prob=05_test_edot -fft -mpi -hdf5 -h5double > configure.log

# Print configuration to stdout
cat configure.log

# Compile
make clean
make -j 16
