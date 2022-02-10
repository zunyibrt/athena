#!/bin/bash

#Load modules
module --force purge
module load modules/1.59-20220201
module load intel-oneapi-compilers/2022.0.1 intel-oneapi-mkl/2022.0.1 intel-oneapi-mpi/2021.5.0 fftw/3.3.10-mpi hdf5/1.12.1-mpi
export MPICH_CXX=icpc

# Run configuration script and save to logfile
./configure.py --prob=clouds_SN -fft -mpi --cxx icpc -hdf5 -h5double > configure.log 

# Print configuration to stdout
cat configure.log

# Compile
make clean
make -j 16
