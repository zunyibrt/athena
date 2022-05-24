#!/bin/bash

#Load modules
module --force purge
module load modules/1.59-20220201
module load intel-oneapi-compilers intel-oneapi-mkl intel-mpi/2017.4.196 fftw/3.3.10-mpi hdf5/1.12.1-mpi
export MPICH_CXX=icpc

# Run configuration script and save to logfile
#./configure.py --prob=07_final --cxx=icpc -b --flux hlld --nscalars=4 -fft -mpi -hdf5 -h5double > configure.log
./configure.py --prob=07_final -b --flux=hlld --nscalars=4 -fft -mpi -hdf5 -h5double > configure.log # Fast compile for testing

# Print configuration to stdout
cat configure.log

# Compile
make clean
make -j 16
