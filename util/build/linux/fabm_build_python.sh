#!/bin/bash

if [ "$1" != "" ]; then
   action=$1
   echo "doing a "$action
else
   action=install
fi

# default Fortran compiler is gfortran - overide by setting compuiler like:
# export compiler=ifort
compiler=${compiler:=gfortran}

host=python

np=-j8
cd $host/$compiler
make $np $action
cd ../..
