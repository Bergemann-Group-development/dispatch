#!/bin/csh

set dir = $0:h
set dir = $dir:h

# Where to find DISPATCH
if (! $?DISPATCH_TOP) then
  if      (-d $dir/config) then
    setenv DISPATCH_TOP $dir
  else if (-d ~/codes/DISPATCH2) then
    setenv DISPATCH_TOP ~/codes/DISPATCH2
  else if (-d ~/codes/dispatch2) then
    setenv DISPATCH_TOP ~/codes/dispatch2
  else if (-d ~/codes/DISPATCH) then
    setenv DISPATCH_TOP ~/codes/DISPATCH
  else if (-d ~/codes/dispatch) then
    setenv DISPATCH_TOP ~/codes/dispatch
  else
    echo "DISPATCH_TOP must be set in the environment"
  endif
endif
echo "DISPATCH_TOP = $DISPATCH_TOP"

# Use the host.pl script to propperly set the HOST variable
set opath = ( $path )
set path = ( ~/bin $DISPATCH_TOP/utilities/scripts $path )
# Set up HOST dependent path
setenv HOST `host.pl`  
if ("$HOST" == "") then
  echo "no HOST env var defined"
  set path = ( ~/bin $DISPATCH_TOP/config/host/$HOST $DISPATCH_TOP/utilities/scripts $path )
else
  echo "HOST = $HOST"
endif
if (-e $DISPATCH_TOP/config/host/$HOST/env.csh) then
  source $DISPATCH_TOP/config/host/$HOST/env.csh
endif
# Use a SCHEDULER variable, if existing
if ($?SCHEDULER) then
  set path = ( ~/bin $DISPATCH_TOP/config/host/$HOST $DISPATCH_TOP/utilities/scripts/$SCHEDULER $DISPATCH_TOP/utilities/scripts $opath )
  echo "SCHEDULER = $SCHEDULER"
endif
