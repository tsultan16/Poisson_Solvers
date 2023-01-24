#!/bin/sh
gfortran -c input.f90
gfortran -O2 input.o SOR.f90  -o SOR -ffpe-trap=invalid,zero,overflow
./SOR

