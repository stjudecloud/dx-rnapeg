#!/bin/bash
# Submits a job array where each job will be taken from a single line of a
# script.
#
# You may also pass --dry-run before the first parameter to simply print out the
# bsub command
#
# You may use -J to specify the job name and/or the maximum number of jobs
# running at once.  If you do specify -J, then this script will attempt to
# parse out and use the job name and/or concurrently running task limit.
#
# Additionally, there are defaults filled in for -o, -e, -M, and 
# -R rusage[mem=?].  If you specify -o and/or -e, then the defaults are not
# used.  If you specify -M without -R rusage[mem=?], then the memory reservation
# (-R rusage[mem=?]) defaults to being the same as -M.  To see defaults,
# run with --dry-run and a commands file, and no optional parameters.
#
# The --mem custom parameter is a platform-neutral alias for -M.
#
# If AFC_DEFAULT_QUEUE is specified in the environment, and no queue is
# specified in your command, then -q $AFC_DEFAULT_QUEUE is included
#
# If AFC_ARGS is set in the environment, then those arguments are passed to
# bsub
#
# $1 = file with one task on each line
# $2,... = other bsub parameters (can include -J; see above)

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
SCRIPT=$1
if [ "$SCRIPT" == "--dry-run" ]; then DRY=$1; shift; SCRIPT=$1; fi
shift
# Do not store the other args in ARGS, as it doesn't work right with quoting
#ARGS="$@"

# Analyze script file
numlines=`cat $SCRIPT | wc -l`
firstline=`head -n 1 $SCRIPT`
if [ "${firstline:0:1}" == "#" ]
then
  start=2
  firstline=`head -n 2 $SCRIPT | tail -n 1`
else
  start=1
fi

# Some defaults
#DEFAULT_MEM=3800
DEFAULT_MEM=
DEFAULT_OEDIR=.

# Look for settings for -P, -J, -o, -e, as these will be defaulted in if not
# specified; also if the environment variable for -P is set, then we will not
# put in the default -P
COMING=_
name=
limit=
hasoe=
oedir=
oetpl=
hasq=
hasmemres=
hasmemlim=
memlim=
res=
resource=
bsub_args=
#for param in $ARGS
while [ "$#" -gt 0 ]
do
  param="$1"
  shift
  
  # Explicitly setting these to their defaults is the signal to read from next
  # parameter.  This is an easy way to handle the case where the switch comes
  # at the end of the parameter list
  if [ "$memlim" == "$COMING" ]; then memlim=$param; fi
  if [ "$oedir" == "$COMING" ]
  then
    oedir=$param
    if [ ! -d "$oedir" ]; then echo "Log dir does not exist: $oedir" >&2; exit 1; fi
    continue
  fi
  if [ "$oetpl" == "$COMING" ]; then oetpl=$param; continue; fi
  if [ "$res" == "$COMING" ]
  then
    res=
    for par in $param
    do 
      resource="$resource -R '$par'"
    done
    if echo "$param" | grep -q mem=
    then hasmemres=1
    fi
    continue
  fi
  
  # Handle -J as a special case, as we always must set/modify -J
  if [ "$name" == "$COMING" ]
  then
    read name limit < <(echo $param | sed 's/\[.*\]//' | awk -F % '{ if($1 == "") $1 = "_"; print $1 " " $2 }')
    continue
  elif [ "$param" == "-J" ]
  then
    name=$COMING
    continue
  fi
  
  # See if -o or -e was set
  if [ "$param" == "-o" ]
  then hasoe=1
  elif [ "$param" == "-e" ]
  then hasoe=1
  fi
  
  # See if -q was set
  if [ "$param" == "-q" ]
  then hasq=1
  fi
  
  # Make --mem an alias for -M
  if [ "$param" == "--mem" ]
  then param="-M"
  fi
  
  # See if -M was set
  if [ "$param" == "-M" ]
  then
    hasmemlim=1
    # Signal that mem lim should be read from next param
    memlim=$COMING
  fi
  
  # See if -R was set
  if [ "$param" == "-R" ]
  then res=$COMING; continue
  fi
  
  # See if custom --log-tpl was set, and if so signal that it should be read
  # from the next param
  if [ "$param" == "--log-dir" ]
  then oedir=$COMING; continue
  elif [ "$param" == "--log-tpl" ]
  then oetpl=$COMING; continue
  fi
  
  # Add to args
  if echo "$param" | grep -qE '[^A-Za-z0-9_-]'
  then bsub_args="$bsub_args '$param'"
  else bsub_args="$bsub_args $param"
  fi
done
# Make sure there wasn't a dangling specification of a tracked option
if [ "$name" == "$COMING" -o "$memlim" == "$COMING" -o "$oedir" == "$COMING" -o "$oetpl" == "$COMING" -o "$res" == "$COMING" ]
then echo "Dangling option with no value" >&2; exit 1
fi
# Get a default name from the command file, if name was not set in args
if [ "$name" == "" ]
then
  name=`echo "$firstline" | awk '{ print $1 }'`
  if [ "$name" == "java.sh" ]
  then
    name=`echo "$firstline" | grep -woE '[^ ]*\.jar' | sed 's/.jar$//'`
    if [ "$name" != "" ]
    then
      name=`basename $name`
    else
      name=`echo "$firstline" | cut -c 8- | sed -r 's/( [^- ][^- ]*).*$/\1/' | sed 's/.*[.]//'`
      if [ "$name" == "" ]; then name="java.sh"; fi
    fi
  else
    name=`basename $name`
  fi
fi

# Add -o and -e if not set
if [ $hasoe ]
then :
else
  if [ $oetpl ]
  then :
  else
    if [ $oedir ]
    then :
    else oedir=$DEFAULT_OEDIR
    fi
    oetpl=$oedir/%N.%J.%I.%W.txt
  fi
  o=`echo $oetpl | sed -e "s/%N/$name/g" -e "s/%W/out/g"`
  e=`echo $oetpl | sed -e "s/%N/$name/g" -e "s/%W/err/g"`
  bsub_args="$bsub_args -o $o -e $e"
fi

# Add -q if not set and AFC_DEFAULT_QUEUE is specified
if [ $hasq ]
then :
elif [ $AFC_DEFAULT_QUEUE ]
then bsub_args="$bsub_args -q $AFC_DEFAULT_QUEUE"
fi

# Add -M if not set
if [ $hasmemlim ]
then :
elif [ "$DEFAULT_MEM" != "" ]
then
  memlim=$DEFAULT_MEM
  bsub_args="$bsub_args -M $memlim"
fi

# Add -R rusage[mem=?] if not set
if [ $hasmemres ]
then :
elif [ "$memlim" != "" ]
then
  resource="$resource -R 'rusage[mem=$memlim]'"
fi

# Add -R
bsub_args="$bsub_args $resource"

# If the command file was empty, then force dry run
if [ "$numlines" -le 0 ]
then echo "Commands file empty; falling back to dry run"; DRY=dry
fi

# Do the bsub(s)
from=$start
if [ "$AFC_SPLIT_CUTOFF" ]
then max=$AFC_SPLIT_CUTOFF
else max=$numlines
fi
while [ $from -le $numlines ]
do
  # Determine high index of range
  let to=from+max-1
  if [ $to -gt $numlines ]; then to=$numlines; fi
  
  # Add the -J param
  j_arg="-J $name[$from-$to]"
  if [ $limit ]; then j_arg="$j_arg%$limit"; fi
  
  # Echo or run command
  cmd="bsub $AFC_ARGS $j_arg $bsub_args lsf_array_for_cmdfile.sh $SCRIPT"
  if [ "$DRY" == "" ]
  then
    echo "Submitting: $cmd"
    set -o pipefail
    eval "$cmd" | tee >( grep 'is submitted to ' | grep -oE '[0-9]*' | head -n 1 > ~/.afc_last_jobid )
    exitcode=$?
  else
    echo "Command: $cmd"
    exitcode=$?
  fi
  
  # If an error was encountered, then break
  if [ "$exitcode" != 0 ]
  then
    if [ "$from" != "$start" ]
    then echo; echo "ALERT: Some jobs were submitted before there was a failure!" >&2
    fi
    break
  fi
  
  # Update low index for next iteration
  let from=to+1
done

# Print timestamp
echo "["`date`"]"

exit $exitcode
