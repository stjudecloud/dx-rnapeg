#!/bin/bash
# Sets up a single notification for an array job.
#
# Submit the array job first (using -o and -e but NOT -N); be sure to
# get the job id.
#
# If you specify 'auto' as the jobid, then the job ID is read from
# ~/.afc_last_jobid, and the message will be started with `pwd` and the contents
# of last_step.txt in the current directory, if it exists
#
# $1 = job id, or 'auto' (see above)
# $2,... = any other message you want to include

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
ID=$1
shift
MORE="$@"

# Handle auto
if [ "$ID" == "auto" ]
then
  if [ ! -s ~/.afc_last_jobid ]
  then echo ~/.afc_last_jobid does not exist, or is empty >&1; exit 1
  fi
  ID=`cat ~/.afc_last_jobid`
  automsg="`pwd` `head -n 1 last_step.txt 2>/dev/null`"
fi

# Do the submission
bsub -M 1000 -P bnotify -w "ended($ID)" -J NOTIFY_$ID echo JOB ENDED: $ID $automsg $MORE
