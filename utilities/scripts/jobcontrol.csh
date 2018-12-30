#!/bin/csh

# job control. Syntax: jobcontrol.csh script [ arguments ]  
#
# where 'script' is an executable script
#
# Touching restart.flag makes the job retry again after
# stopping (forced with 'touch stop.flag' for stagger-code).
#
# To stop the job from looping with failure, remove the
# restart.flag manually.  The job waits a maximum of retrymax
# seconds.

#set echo							# for debugging

# Where to find DISPATCH
if (! $?DISPATCH_TOP) then
  if      (-d ~/codes/DISPATCH2) then
    setenv DISPATCH_TOP ~/codes/DISPATCH2
  else if (-d ~/codes/dispatch2) then
    setenv DISPATCH_TOP ~/codes/dispatch2
  else if (-d ~/codes/DISPATCH) then
    setenv DISPATCH_TOP ~/codes/DISPATCH
  else if (-d ~/codes/dispatch) then
    setenv DISPATCH_TOP ~/codes/dispatch
  else
    echo "DISPATCH_TOP must be set in the environment"
    exit
  endif
endif

set path = ( $DISPATCH_TOP/utilities/scripts $DISPATCH_TOP/utilities/scripts/slurm $path )

foreach f ( jobrunning jobkill jobid jobmail jobhelper )	# check that commands are in path
  which $f >/dev/null 
  if ($status) then
    echo "command $f is not in path; adjust your ~/.cshrc or ~/.tcshrc"
    echo "csh path = $path"
    exit 1
  else
    echo "using `which $f`"
  endif
end

if (-e redir.flag) then						# redir task
  set d = `cat redir.flag`					# read
  \rm redir.flag						# remove
  cd $d								# go there
endif

echo "master: `hostname`"					# show master node name

echo "cmd = $*"

#---------------------------------------------------------------
@ m = 0								# mail times
@ n = 0 							# attempts
@ s = 5								# check time
@ j = 240							# grace period for stopping
if (-r jobstop.flag) @ j = `cat jobstop.flag` 			# non-default time to stop

@ wt = 0
start:

set jid = `jobid`
set oid = "nojob"						# default value
if (-r xunning.flag) then					# another jobs is running here
  if ($wt == 0) then					        # first time
    if (-r jobstop.flag) set j = `cat jobstop.flag`		# non-default time to stop
    echo "touching dump.flag, waiting max $j seconds"		# info
    touch takeover.flag                 			# mark takeover as ongoing
    touch dump.flag 						# ask oid to dump & stop
  else if ($wt > $j) then					# we tried the nice way
    echo "jobkill $oid after $wt seconds"			# info
    set oid = `cat xunning.flag`				# the other jobid
    if ($oid == $jid) then                                      # check suicide
      echo "refusing to kill myself"                            # don't
    else                                                        # bondafide
      jobrunning $oid						# check if it runs
      if ($status) then						# apparently not
        echo "$oid is apparently not running, removing running.flag"# info
        \rm running.flag					# not running
      else							# running
        echo "$oid is apparently running, issuing jobkill"	# info
        jobkill $oid						# do it the hard way
      endif                                                     # end of kill attempt
    endif                                                       # end of wait
    goto cont
  endif
  sleep $s							# wait briefly
  @ wt += $s							# accumulated waiting time
  echo "$wt / $j"                                               # info
  goto start							# try again
endif

cont:
if ($wt > 0) echo "waited $wt seconds"				# info

#jobid > running.flag						# set our running flag
\rm -f dump.flag stop.flag takeover.flag			# remove the flags we set

#---------------------------------------------------------------
@ t =  0							# retry time
@ m =  0							# mail times
@ s = 120							# process restar
@ retrymax = 18000						# 5h

( cd data/$2; jobhelper >& jobhelper.log </dev/null & )		# start jobhelper
if (-e jobhelper.flag) \rm jobhelper.flag			# remove any initial flag

restart:

@ w = `date +%s`						# UNIX seconds
echo "$*"
$*								# user script, local or in ~/bin
@ w = `date +%s` - $w						# run time in sec

if ($w < 0) then						# should not happen
  @ w = $s							# wait nominal
else if ($s > $w) then						# short run time
  @ w = $s - $w							# reduced wait
else								# long run time
  @ w = 0							# zero wait
endif

if (-e restart.flag) then					# should we attempt to restart?
  #\rm -f restart.flag						# repear or not?
  @ m++
  if ($m < 10 || $m == 50 ) then				# first 9 times, and 50th time
    echo "$jid restarted $m times" | jobmail $user		# mail user on restart
  endif
  if (-e retrymax.flag) then					# allow change
    @ retrymax = `cat retrymax.flag`
    \rm retrymax.flag
  endif
  if (-e data/$2/jobhelper.flag) then				# new jobhelper flag?
    killall jobhelper						# kill any of ours
    ( cd data/$2; jobhelper >& jobhelper.log </dev/null & )	# start jobhelper
    \rm -f data/$2/jobhelper.flag				# clear flag
  endif
  if (-e wait.flag) then					# allow change
    @ s = `cat wait.flag`					# new wait time
    @ w = $s							# immediate effect
    \rm wait.flag						# remove
    @ wt = $w
  endif
  if (! -e takeover.flag) then					# unless another job takes over
    if ($w > 0) then						# need to wait?
      @ wt = 0							# total wait time
      while ($wt < $w)
        if (-e wait.flag) then					# allow change
          @ s = `cat wait.flag`					# new wait time
          @ w = $s						# immediate effect
          \rm wait.flag						# remove
	  @ wt = $w
        endif
        ( sleep 5 )						# short wait, incrementally
        @ wt += 5						# accumulate short wait
      end
    endif
    if (-e cd.flag) then					# cd.flag during loop (not takeover)
      set d = `cat cd.flag`					# read
      \rm cd.flag						# remove
      cd $d							# go there
      @ t = 0							# reset retry time
      touch restart.flag					# inherit
      goto restart						# restart
    endif
    @ t += $s							# accumulated retry time
    echo "$t / $retrymax"                                       # info
    if ( $t < $retrymax) goto restart				# don't wait too long
  endif
endif

\rm -f running.flag						# clean exit

if (! -e takeover.flag) then					# if we are not being taken over
  if (-e next.flag) then					# next task (no restart.flag)
    set d = `cat next.flag`					# read
    # \rm next.flag						# remove
    cd $d							# go there
    @ t = 0							# reset retry time
    goto restart						# restart
  endif
endif

killall jobhelper
