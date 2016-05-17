#!/bin/sh

# if not set use the suggested source code installation directories
FABM_BASE=${FABM_BASE:=~/FABM/code}
GOTM_BASE=${GOTM_BASE:=~/GOTM/code}

# default Fortran compiler is gfortran - overide by setting compuiler like:
# export compiler=ifort
compiler=${compiler:=gfortran}

host=python
# ready to configure
mkdir -p $host/$compiler
cd $host/$compiler
cmake $FABM_BASE/src/drivers/python \
      -DFABM_EMBED_VERSION=on \
      -DCMAKE_Fortran_COMPILER=$compiler
cd ../..
