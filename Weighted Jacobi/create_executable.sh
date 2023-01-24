#!/bin/sh
gfortran -c input.f90
gfortran -O2 input.o WJacobi.f90  -o WJacobi -ffpe-trap=invalid,zero,overflow
./WJacobi

