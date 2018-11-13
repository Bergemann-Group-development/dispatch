#!/bin/bash

# job control. Syntax: jobcontrol.sh script [ arguments ]  
#
# where 'script' is an executable script
#
# Touching restart.flag makes the job retry again after
# stopping (forced with 'touch stop.flag' for stagger-code).
#
# To stop the job from looping with failure, remove the
# restart.flag manually.  The job waits a maximum of retrymax
# seconds.

#set -x

export DISPATCH_TOP=$HOME/codes/dispatch2

export PATH="${DISPATCH_TOP}/config/host/${HOST}:${DISPATCH_TOP}/utilities/scripts:${PATH}"

for f in jobrunning jobkill jobid jobmail jobhelper; do
  if { ! which $f >/dev/null 2>&1; }; then
    echo "============================================================"
    echo "command $f is not in path; adjust your ~/.bashrc"
    echo "============================================================"
    echo "PATH = $PATH"
    exit 1
  fi
done

if [ -e redir.flag ]; then
  d=`cat redir.flag`					        # read
  \rm redir.flag						# remove
  cd $d								# go there
fi

echo "master: `hostname`"					# show master node name
echo "cmd = $*"

#---------------------------------------------------------------
m=0								# mail times
n=0 							        # attempts
s=5								# check time
j=240							        # grace period for stopping
if [ -e jobstop.flag ]; then j=`cat jobstop.flag`; fi           # non-default grace period

jid=`jobid`
oid="nojob"     						# default value
wt=0
if [ -e running.flag ]; then
  if [ -e jobstop.flag ]; then
    j=`cat jobstop.flag`
  fi
  echo "touching dump.flag, waiting max $j seconds"		# info
  touch takeover.flag                 			        # mark takeover as ongoing
  touch dump.flag 						# ask oid to dump & stop
  while { jobrunning $oid; }; do
    echo "$oid is apparently running"				# info
    if [ $wt -gt $j ]; then						# we tried the nice way
      echo "jobkill $oid after $wt seconds"			# info
      if [ $oid == $jid ]; then
        echo "refusing to kill myself"
      else
        jobkill $oid						# do it the hard way
      fi
      /bin/rm running.flag
      break
    fi
    sleep $s							# wait briefly
    let wt=wt+s
  done
fi
if [ $wt -gt 0 ]; then echo "waited $wt seconds"; fi		# info

jobid > running.flag						#ourrunningflag
\rm -f dump.flag stop.flag takeover.fag				#removetheflagswe

#---------------------------------------------------------------
t=0		        					# retry time
m=0		        					# mail times
s=120		        					# process restar
retrymax=18000	        					# 5h

jobhelper >jobhelper.log 2>&1 </dev/null &			# start jobhelper
if [ -e jobhelper.flag ]; then \rm jobhelper.flag; fi		# remove any initial flag

let s=1                                                         # signal number
while [ $s -lt 65 ]; do                                         # for all signals
  trap '' $s                                                    # set insensitive
  let s=s+1                                                     # next signal
done

w=`date +%s`	        					# UNIX seconds
$*								# user script, local or in ~/bin
w=$((`date +%s` - $w))						# run time in sec

let s=1                                                         # signal number
while [ $s -lt 65 ]; do                                         # for all signals
  trap - $s                                                     # set sensitive
  let s=s+1                                                     # next signal
done

if [ $w -lt 0 ]; then						# should not happen
  let w=s							# wait nominal
else
  if [ $s -gt $w ]; then					# short run time
    let w=s-w							# reduced wait
  else								# long run time
    w=0   							# zero wait
  fi
fi

while [[ -e restart.flag && $t -lt $retrymax ]]; do
  let m+=1
  if [[ $m -lt 10 || $m -eq 50 ]]; then				# first 9 times, and 50th time
    echo "$jid restarted $m times" | jobmail $user		# mail user on restart
  fi
  if [ -e retrymax.flag ]; then					# allow change
    retrymax=`cat retrymax.flag`
    \rm retrymax.flag
  fi
  if [ -e jobhelper.flag ]; then				# new jobhelper flag?
    killall jobhelper						# kill any of ours
    jobhelper >jobhelper.log 2>&1 </dev/null &a			# start a new one
    \rm jobhelper.flag						# clear flag
  fi
  if [ -e wait.flag ]; then					# allow change
    s=`cat wait.flag`   					# new wait time
    let w=s							# immediate effect
    #\rm wait.flag						# remove
    let wt=w
  fi
  if [ ! -e takeover.flag ]; then				# unless another job takes over
    if [ $w -gt 0 ]; then					# need to wait?
      wt=0							# total wait time
      while [ $wt -lt $w ]; do
        if [ -e wait.flag ]; then				# allow change
          s=`cat wait.flag`					# new wait time
          let w=s						# immediate effect
          #\rm wait.flag					# remove
	  let wt=w
        fi
        ( sleep 5 )						# short wait, incrementally
        let wt+=5						#accumulateshortwait
      done
    fi
    if [ -e cd.flag ]; then					# cd.flag during loop (not takeover)
      d=`cat cd.flag`   					# read
      \rm cd.flag						# remove
      cd $d							# go there
      t=0							#reretrytime
      touch restart.flag					# inherit
    fi
    let t+=s							# accumulated retry time
    if [ $t -gt $retrymax ]; then exit; fi			# don't wait too long
  fi
done

\rm -f running.flag						# clean exit

killall jobhelper                                               # stop background jobhelper
