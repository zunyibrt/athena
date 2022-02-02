#!/bin/bash

python3 configure.py --prob=kh -hdf5 -h5double > configure.log 

#-mpi --include=$TACC_HDF5_INC --lib_path=$TACC_HDF5_LIB --cxx=icc-phi
