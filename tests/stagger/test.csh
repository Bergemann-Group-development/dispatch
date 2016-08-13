#!/bin/csh

cd $PBS_O_WORKDIR
setenv KMP_AFFINITY "granularity=core,scatter"

make clean
make COMPILER=ifort OPTS=optimized info
make COMPILER=ifort OPTS=optimized
./test.x

make clean
make COMPILER=ifort OPTS=optimized OPT=-O2 HW=-xHost info
make COMPILER=ifort OPTS=optimized OPT=-O2 HW=-xHost
./test.x

make clean
make COMPILER=ifort OPTS=optimized OPT=-Ofast HW=-xHost info
make COMPILER=ifort OPTS=optimized OPT=-Ofast HW=-xHost
./test.x
