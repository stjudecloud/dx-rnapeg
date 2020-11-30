#!/bin/bash
# Submits an array job where each task will be taken from a single line of a
# script.
#
# You may also pass --dry-run before the first parameter to simply print out the
# qsub command
#
# $1 = file with one task on each line
# $2,... = other qsub parameters (exclude -t and -b; can exclude -N)

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
SCRIPT=$1
if [ "$SCRIPT" == "--dry-run" ]; then DRY=$1; shift; SCRIPT=$1; fi
shift
ARGS=$*

# Analyze script file and parameters
numlines=`cat $SCRIPT | wc -l`
firstline=`head -n 1 $SCRIPT`
if [ "${firstline:0:1}" == "#" ]
then
  start=2
  firstline=`head -n 2 $SCRIPT | tail -n 1`
else
  start=1
fi
hasname=no
for param in $*
do
  if [ "$param" == "-N" ]
  then hasname=yes ; break
  fi
done
if [ "$hasname" == "no" ]
then
  name=`echo "$firstline" | awk '{ print $1 }'`
  if [ "$name" == "java.sh" ]
  then
    name=`echo "$firstline" | grep -oE '[^ ]*.jar' | sed 's/.jar$//'`
    if [ "$name" == "" ]
    then
      name=`echo "$firstline" | cut -c 8- | sed -r 's/( [^- ][^- ]*).*$/\1/' | sed 's/.*[.]//'`
      if [ "$name" == "" ]; then name="java.sh"; fi
    fi
  fi
  nameargs="-N $name"
fi

# Do the qsub
cmd="qsub $nameargs -b y -t $start:$numlines $ARGS sge_array_for_cmdfile.sh $SCRIPT"
echo "Submitting: $cmd"
if [ "$DRY" == "" ]
then $cmd
fi

# Warn if no args
if [ "$ARGS" == "" ]; then echo "WARNING: NO QUEUE SPECIFIED"; fi