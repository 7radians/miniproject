#!/bin/bash

# Modified from Workshop 7 build script

# Fortran files to compile
programfiles="command_line.f90 create_axis.f90 particle_mover.f90 write_netcdf.f90 gauss.f90"

# Name of compiled file
outfile="miniproject"

# Compiler and compiler flags
fc=gfortran
fflags=`nf-config --fflags`
flibs=`nf-config --flibs`

# Compile line
$fc -g -std=f2008 $fflags $programfiles $flibs -o $outfile -Wall
if [ $? -ne 0 ]
then
  echo Compilation error, terminating.
  return 1
else
  echo Compilation successful.
fi

# Run fortran with given input data
time ./$outfile nx=5 ny=5 problem=null
if [ $? -ne 0 ]; then
  echo Error generating data, terminating.
  return 1
else
  echo Data successfully generated.
fi

# Run visualisation code
time python3 readfile.py
if [ $? -ne 0 ]; then
  echo Error visualising data, terminating.
  return 1
else
  echo Data visualisations successfully generated.
fi