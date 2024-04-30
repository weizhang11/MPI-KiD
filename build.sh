#!/bin/bash

make clean

ulimit -s unlimited
ulimit -v unlimited

#make COMPILER=nvfortran DEF_CASIM=1 CASE=SC_2D all
#make -j 20 VENDOR=nvidia COMPILER=ftn PCOMPILER=ftn DEF_CASIM=1 DEF_DPREC=true CASE=SC_2D all
make -j 20 VENDOR=cray COMPILER=ftn PCOMPILER=ftn DEF_CASIM=1 DEF_DPREC=true CASE=CU_2D all


#make CASE=1D all
