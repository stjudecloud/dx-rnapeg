#/bin/bash
#
# This script periodically checks whether the given LSF job is complete. 
# If it's not, it sleeps; if it's complete, it exits. Following are the exit 
# codes this script uses:
#   0: all jobs completed successfully
#   1: any error not covered by other non-zero exit codes mentioned below
#   2: at least one job in the job array failed
#   3: job status of the given job ID could not be found 
#
# Parameters:
# $1 = LSF job ID
# $2 = Sleep interval in seconds (Ex: 10s)
# $3 = (optional) file to write bjobs output to if there is a failure

# Show usage information if no parameters were sent
if [ $# -lt 2 ]; then about.sh $0; exit 1; fi

# Get parameters
lsf_job_id=$1
sleep_interval=$2
failed_bjobs_out=$3

# Exit if job status could not be found.
last_status=`bjobs $lsf_job_id 2>/dev/null`
if [ -z "$last_status" ]; then exit 3; fi

# bjobs output file for debugging purposes
bjobs_out=`mktemp`

# Loop until all jobs are done. Sleep between querying bjobs.
while :
do
  # Check job status
  last_status=`bjobs $lsf_job_id 2>/dev/null | tee $bjobs_out | tail -n +2 | awk '{print $3}' | sort | uniq | tail -n 1`
  
  # All jobs are completed successfully
  if [ "$last_status" == "DONE" ]
  then
    rm $bjobs_out
    exit 0
  # At least one job has failed
  elif [ "$last_status" == "EXIT" ]
  then
    if [ "$failed_bjobs_out" == "" ]
    then rm $bjobs_out
    else mv $bjobs_out $failed_bjobs_out
    fi
    exit 2
  fi
  
  # The jobs haven't completed, so sleep for a while
  sleep $sleep_interval
done
