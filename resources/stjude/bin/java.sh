#!/bin/bash
# Calls java, including some default arguments if not specified.
#
# Defaults:
# If you do not specify -Xmx, then -Xmx2g is used
# If you do not specify -XX:ParallelGCThreads, then -XX:ParallelGCThreads=4 is used
#
# $1... = all arguments are passed to java command

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
ARGS="$@"

# Default options
hasxmx=
hasgct=
for arg in $ARGS
do
  if [ "${arg:0:4}" == "-Xmx" ]; then hasxmx=1
  elif [ "${arg:0:22}" == "-XX:ParallelGCThreads=" ]; then hasgct=1
  fi
done
if [ $hasxmx ]; then : ; else ARGS="-Xmx2g $ARGS"; fi
#if [ $hasgct ]; then : ; else ARGS="-XX:ParallelGCThreads=1 $ARGS"; fi

# DEBUG
# AMF - 27JUN12
# MCR - 23JUL12
#MEM_ARGS="-Xmx4g -Xms4g -Xmn2560m -Xss160k"
#DEBUG_ARGS="-XX:-PrintCommandLineFlags -XX:+PrintGCDetails -XX:+PrintGCTimeStamps -XX:+TraceClassLoading"
JAVA_TMPDIR=`mktemp -d`
OPT_ARGS="-Djava.io.tmpdir=$JAVA_TMPDIR -XX:ParallelGCThreads=3 -XX:+UseConcMarkSweepGC -XX:+UseParNewGC -XX:SurvivorRatio=8 -XX:TargetSurvivorRatio=90 -XX:MaxTenuringThreshold=15"
#if [ "`whoami`" == "pgupta" ]; then OPT_ARGS="-XX:ParallelGCThreads=1 -XX:+UseConcMarkSweepGC"; fi
ARGS="${OPT_ARGS} ${ARGS}"
#ARGS="${MEM_ARGS} ${DEBUG_ARGS} ${OPT_ARGS} ${ARGS}"

# Remote debugging
if [ "$JAVA_REMOTE_DEBUG" != "" -a "$JAVA_REMOTE_DEBUG_ENABLED" != "" ]
then ARGS="-Xdebug -Xrunjdwp:transport=dt_socket,address=$JAVA_REMOTE_DEBUG,suspend=y $ARGS"
fi

# Print out command
if [ $JAVA_VERBOSE ]; then echo java $ARGS; fi

# Call java
if [ $JAVA_POLLING_PMAP ]
then
  java $ARGS &
  PID=$!
  polling_pmap.sh $PID 5 120 >&2
  wait $PID
else
  java $ARGS
fi
