#!/bin/bash
# Runs pmap on a process every n seconds
# Exits once that process no longer exists
# Output is written to stderr
#
# $1 = process id
# $2 = polling interval in seconds (default 60)
# $3 = reporting interval in number of polls (default 1 for every poll)

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
PID=$1
POLL_INTERVAL=$2
RPT_POLLS=$3
if [ "$POLL_INTERVAL" == "" ]; then POLL_INTERVAL=60; fi
if [ "$RPT_POLLS" == "" ]; then RPT_POLLS=1; fi

# Get 2 temporary files and a variable to help with swapping pointers to the
# two files
cur=`mktemp`
last=`mktemp`
swapper=

# Register trap
trap "echo 'Caught signal; last two polls in any order:'; cat $cur; echo; cat $last; exit" TERM

# Do poll loop
polls=0
totalthresh=1000000
pollthresh=0
reported=
while pmap -x $PID > $cur
do
  # Add footer to report
  echo "host=`hostname`, PID=$PID INTERVAL=$POLL_INTERVAL poll #$polls, `date`" >> $cur
  
  # See if we need to report
  total=`grep total $cur | tr -d [:alpha:][:space:]`
  if [ $polls -ge $pollthresh -o $total -ge $totalthresh ]
  then
    cat $cur
    echo
    reported=1
    let "pollthresh = polls + RPT_POLLS"
    let "totalthresh = total + 200000"
  else
    reported=
  fi
  
  # Swap cur and last (to assign last and prepare cur for next run
  swapper=$last
  last=$cur
  cur=$swapper
  # Increment # polls
  let ++polls
  
  # Sleep
  sleep $POLL_INTERVAL
done

# Make sure the last poll's results are reported
if [ ! $reported ]
then
  cat $last
  echo
fi

# Report completion
echo "Process gone"
echo "host=`hostname`, PID=$PID INTERVAL=$POLL_INTERVAL poll #$polls, `date`"
