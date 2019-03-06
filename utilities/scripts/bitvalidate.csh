#!/bin/csh

# Check if the input file exists
set nml = $1/bitvalidate.nml
if (-e $nml) then
  echo "\n$nml OK"
else
  mkdir -p $1 data/$1/bitvalidate
  cp bitvalidate.nml $nml
endif
git add -f $nml

# Check if the validation file exists
set val = data/$1/bitvalidate/rank_00000.val
if (-e $val) then
  echo "\n$val OK"
  make --no-print-directory SOLVER=$1 OPTS=bitvalidate -j
else
  cp -p $nml $nml.bak
  make clean
  make --no-print-directory SOLVER=$1 OPTS=bitvalidate -j
  perl -i -p -e "s{mode='compare'}{mode='write'}" $nml
  setenv OMP_NUM_THREADS 1
  ./dispatch.x $nml
  mv $nml.bak $nml
endif
sleep 3
