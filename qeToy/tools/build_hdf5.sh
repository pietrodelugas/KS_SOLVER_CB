#!/bin/bash
# Builds a serial HDF5 (C with gcc, Fortran with nvfortran of the given HPC SDK) in ~/opt/hdf5-<ver>-nvhpc<sdk>; needed for ham.h5 (cmake -DQE_ENABLE_HDF5=ON -DHDF5_ROOT=<prefix> -DHDF5_USE_STATIC_LIBRARIES=TRUE)
# usage: build_hdf5.sh <hpc_sdk_version e.g. 26.3>  -> installs serial HDF5 (C by gcc, Fortran by nvfortran) in ~/opt/hdf5-1.14.6-nvhpc<ver>
set -e
V=1.14.6; SDK=$1; S=$(dirname "$0"); PREFIX=$HOME/opt/hdf5-$V-nvhpc$SDK
export PATH=${HPCSDK_ROOT:-/opt/nvidia/hpc_sdk}/Linux_$(uname -m)/$SDK/compilers/bin:$PATH
cd "$S"
[ -f hdf5-$V.tar.gz ] || curl -sL -o hdf5-$V.tar.gz https://github.com/HDFGroup/hdf5/releases/download/hdf5_$V/hdf5-$V.tar.gz
rm -rf src_$SDK && mkdir src_$SDK && tar xzf hdf5-$V.tar.gz -C src_$SDK 
cd src_$SDK/hdf5-$V && mkdir b && cd b
CC=gcc FC=nvfortran ../configure --prefix=$PREFIX --enable-fortran --disable-shared --disable-tests --disable-tools > conf.log 2>&1
make -j8 > make.log 2>&1
make install > inst.log 2>&1
echo "installed in $PREFIX"; ls $PREFIX/include | grep -i "hdf5.mod"
