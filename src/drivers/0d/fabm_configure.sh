#!/bin/sh

# if not set use the suggested source code installation directories
FABM_BASE=${FABM_BASE:=~/FABM/code}
GOTM_BASE=${GOTM_BASE:=~/GOTM/code}

# default Fortran compiler is gfortran - overide by setting compuiler like:
# export compiler=ifort
compiler=${compiler:=gfortran}

# NetCDF
# nf-config must be in the path and correpsond to the value of compiler
# try:
# nf-config --all

# Make install prefix configurable
install_prefix=${install_prefix:=~/local/fabm}

host=0d
# ready to configure
mkdir -p $host/$compiler
cd $host/$compiler
cmake $FABM_BASE/src/drivers/0d \
      -DFABM_EMBED_VERSION=on \
      -DGOTM_BASE=$GOTM_BASE/ \
      -DCMAKE_Fortran_COMPILER=$compiler \
      -DCMAKE_INSTALL_PREFIX=$install_prefix/$compiler
cd ../..
