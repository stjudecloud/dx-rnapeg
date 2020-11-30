#!/usr/bin/env bash
#
# This script is a wrapper script that wraps a lower-level script 
# and additionally contains code to print helpful messages and  
# send notifications. 
#
# The wrapped lower-level script periodically checks whether the given 
# LSF job is complete. If it's not, it sleeps; if it's complete, it 
# exits. Following are the exit codes this wrapped lower-level script uses:
#   0: all jobs completed successfully
#   1: any error not covered by other non-zero exit codes mentioned below
#   2: at least one job in the job array failed
#   3: job status of the given job ID could not be found 
#
# USAGE:
# wait_until_job_completes.sh -j STR -s STR -d DIR [options]
#
# PARAMETERS:
# -j STR    List of comma-separated LSF job IDs
#           (Ex: 8765432 or 8765432,8765433,8765434).
# -s STR    Step number (Ex: 05).
# -d DIR    Path to the run dir.
#
# OPTIONS:
# -i STR  (Optional) Sleep interval in seconds (Default: 10s).
# -n STR  (Optional) List of comma-separated email addresses of people 
#         who should be notified of job completion or failure 
#         (Ex: john.doe@stjude.org,jane.smith@stjude.org).
#         The value of the environment variable USER_EMAIL is appended 
#         to this list.     
# -h STR  (Optional) Show usage.

# Show usage information if no parameters were sent
if [ "$#" -lt 6 ]; then about.sh $0; exit 1; fi

# Notify people of job failure
# $1 = Recipient email addresses, comma-separated 
# $2 = Subject
# $3 = From
# $4 = Message
send_email() {
  local from=$1
  local to=$2
  local subject=$3
  local message=$4
  echo $message | mail -s "$subject" -r $from $to 
}

# Wraps send_email function
# $1 = Subject
# $2 = Message
notify() {
  if [ "$notify_recipients" ]
  then
    local subject=$1
    local message=$2
    send_email $USER "$notify_recipients" "$subject" "$message"
  fi
}

# Get parameters
job_id=
step=
run_dir=
sleep_interval="10s"
notify_recipients=
while getopts "hj:s:d:i:n:" OPTION
do
  case $OPTION in
    h)
      about.sh $0
      exit 1
      ;;
    j)
      job_id=$OPTARG
      ;;
    s)
      step=$OPTARG
      ;;
    d)
      run_dir=$OPTARG
      ;;
    i)
      sleep_interval=$OPTARG
      ;;
    n)
      notify_recipients=$OPTARG
      ;;
  esac
done

# Make sure job ID(s) are specified
if [ "$job_id" == "" ]
then
  echo "Error: Job ID is required, but it has not been specified!"
  exit 1
elif [ "$step" == "" ]
then
  echo "Error: Step is required, but it has not been entered!"
  exit 1
elif [ "$run_dir" == "" ]
then
  echo "Error: Run dir path is required, but it has not been entered!"
  exit 1
fi

# Add the user email in the notification list
if [ "$USER_EMAIL" ]
then
  if [ "$notify_recipients" ]
  then
    delim=","
    notify_recipients=${notify_recipients}${delim}${USER_EMAIL}
    notify_recipients=`deduplicate_list.sh $notify_recipients $delim`
  else
    notify_recipients=${USER_EMAIL}
  fi
fi

# Call the lower-level wrapped script
space_separated_list=`echo $job_id | tr "," " "`
check_lsf_job_status.sh "$space_separated_list" $sleep_interval
exit_code=$?
if [ "$exit_code" -ne 0 ]
then
  if [ "$exit_code" -eq 2 ]
  then
    echo "Error: At least one job in the job array $job_id failed! Exiting..."
    subject="Error in wait_until_job_completes.sh $job_id"
    message="Error: At least one job in the job array $job_id failed!"
    notify "$subject" "$message"
    exit $exit_code
  elif [ "$exit_code" -eq 3 ]
  then
    echo "Error: Job status of at least one of the given job IDs ($job_id) could not be found! Exiting..."
    subject="Error in wait_until_job_completes.sh $job_id"
    message="Error: Job status of at least one of the given job IDs ($job_id) could not be found!"
    notify "$subject" "$message"
    exit $exit_code
  else
    echo "Error: check_lsf_job_status.sh failed! Exiting..."
    subject="Error in wait_until_job_completes.sh $job_id"
    message="Error: wait_until_job_completes.sh failed!"
    notify "$subject" "$message"
    exit $exit_code
  fi
fi

# Go through the stdout files and make sure none of the jobs failed. This 
# is necessary because sometimes it's possible to miss detecting a failed 
# job by just using bjobs.
# Before checking the log files, sleep for a while. This is necessary because
# apparently LSF changes the status in the log file "after" reporting DONE 
# status in bjobs.
sleep $sleep_interval
log_dir=$run_dir/file*-$step*/logs
echo "Checking stdout in $log_dir to see if any jobs failed..."
for i in $space_separated_list;
do
  # Make sure there are log files for this job id in the log dir
  num_log_files=`ls -1 $log_dir/out.${i}.*.txt | wc -l`
  if [ "$num_log_files" -eq 0 ]
  then 
    echo "Error: There are no log files for job id ${i} in log dir ${log_dir}."
    exit 1
  fi
  
  # Check the stdout files for any job failures
  result=`grep -L "Success" $log_dir/*out.${i}.*.txt`
  if [ ! -z "$result" ]
    then
    echo "`date`; Step$step; Error: At least one job failed"
    subject="Error in wait_until_job_completes.sh for job ID ${i}"
    message="Error: At least one job in the job array ${i} failed! Step number is $step and run dir is $run_dir"
    notify "$subject" "$message"
    exit 2
  fi
  done
  echo "...No jobs failed in this step!"

# Print success message
echo "Success! Job(s) $job_id completed successfully!"
