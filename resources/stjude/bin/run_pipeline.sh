#!/bin/bash
# Runs each step of a pipeline.
#
# Parameters:
# -d STR  Pipeline run directory path.
# -t STR  (Optional) Pipeline type (Ex: coverage, conserting, etc.).
# -i STR  (Optional) Sleep interval to be used in periodically checking 
#           job status (Ex: 1m (default), 15m, 3h, etc).
# -s NUM  (Optional) Pipeline step from where pipeline should be started (Ex: 01).
# -m NUM  (Optional) Maximum pipeline step to run.
# -f NUM  (Optional) Force run the given step (Ex: 01). If -s is specified, 
#           the value of -f must not be less than the value of -s.
# -n STR  (Optional) List of comma-separated email addresses of people who should 
#           be notified of job completion or failure 
#           (Ex: john.doe@stjude.org,jane.smith@stjude.org).
#           The value of the environment variable USER_EMAIL is appended to this 
#           list.
# -N      (Optional) Notify for job failures only.
# -l STR  (Optional) Label or name for the pipeline. If this param is not 
#           specified, the run dir basename is used as the label.
# -h      (Optional) Show usage.

# Show usage information if no parameters were sent
if [ "$#" -lt 2 ]; then about.sh $0; exit 1; fi

# Create a tmp file
tmp_file=`mktemp`
log_backup=`mktemp`

trap cleanup_on_signal SIGHUP SIGINT SIGTERM

cleanup_on_signal() {
  sig=$(($? - 128))
  if [ "$lsf_job_id" ]
  then
    log "Attempting to kill LSF job: $lsf_job_id" 
    bkill $lsf_job_id
  fi
  failure "Caught signal $(kill -l $sig). LSF jobs may be running $lsf_job_id"
}

# Cleans up acquired resources
cleanup() {
  # Remove tmp file
  rm $tmp_file
  rm $log_backup
}

# Log something
# Usage:
#   long_message_generator | log -
#   log "short message"
# $1 = message or - (date and step will be added and can be ommitted from message)
log() {
  local log_file=$run_dir/$pipeline_progress_file
  if [ ! -d `dirname $log_file` ]; then log_file=$log_backup; fi
  
  if [ "$1" == "-" ]
  then tee -a $log_file
  else echo "`date`; Step$step; $1" | tee -a $log_file
  fi
  
  if [ "$log_file" != "$log_backup" ]; then cp $log_file $log_backup; fi
}

# Notify people of job completion or failure
# $1 = Recipient email addresses, comma-separated 
# $2 = Subject
# $3 = From
# $4 = Message
send_email() {
  local from="$1"
  local to="$2"
  local subject="$3"
  local message="$4"
  ( echo -e "Run directory: $run_dir\n\n$message\n\nLog:"; cat $log_backup ) | mail -s "$subject" -r $from $to 
}

# Wraps send_email function
# $1 = Subject
# $2 = Message
notify() {
  if [ "$notify_recipients" -a "$do_notify" ]
  then
    local subject="$1"
    local message="$2"
    send_email $USER "$notify_recipients" "$subject" "$message"
  fi
}

# Perform failure logging and notifications
# $1 = additional failure message
notify_failure() {
  log "FAILED $1"
  notify "[rp] Fail $pipeline_type $label $step" "$1"
}

# Perform success logging and notifications
# $1 = additional success message
notify_success() {
  log "COMPLETED $1"
  notify "[rp] Done $pipeline_type $label" "$1"
}

preprocess() {
  external_script $preprocess
}

# $1 = notify message
# $2 = exit code (default 1)
failure() {
  notify_failure "$1"
  cleanup
  external_script $failure
  local code="$2"
  if [ "$code" == "" ]; then code=1; fi
  exit $code
}

completion() {
  if [ ! "$failure_notification_only" ]
  then notify_success
  fi
  cleanup
  external_script $completion
}

external_script() {
  local script=$1
  if [ -e "$script" ] 
  then
    run_script $script
  elif [ -e "$run_dir/$script" ]
  then
    run_script $run_dir/$script 
  fi
}

run_script() {
  local script=$1
  eval $(readlink -f $script) 
  exitcode=$?
  if [ $exitcode == 126 ];
  then
    echo "Trying to recover from permission denied by running as a bash script" >&2
    bash $script
    exitcode=$?
  fi
}

# Get parameters
run_dir=
run_dir_name=
pipeline_type=
sleep_interval=
start_step=
max_step=
force_step=
notify_recipients=
failure_notification_only=
label=
do_notify=1
preprocess=".preprocess.sh"
failure=".failure.sh"
completion=".completion.sh"
while getopts "hd:t:i:s:f:n:Nl:m:p:e:c:q" OPTION
do
  case $OPTION in
    h)
      about.sh $0
      exit 1
      ;;
    d)
      run_dir=$OPTARG
      run_dir_name=`basename $run_dir`
      ;;
    t)
      pipeline_type=$OPTARG
      ;;
    i)
      sleep_interval=$OPTARG
      ;;
    s)
      start_step=$OPTARG
      ;;
    f)
      force_step=$OPTARG
      ;;
    n)
      notify_recipients=$OPTARG
      ;;
    N)
      failure_notification_only="failure_notification_only"
      ;;
    l)
      label=$OPTARG
      ;;
    m)
      max_step=$OPTARG
      ;;
    p)
      preprocess=$OPTARG
      ;;
    e)
      failure=$OPTARG
      ;;
    c)
      completion=$OPTARG
      ;;
    q)
      do_notify=
      ;; 
  esac
done

if [[ -z "$run_dir_name" ]]; then
    echo "You must specify a run directory (-d)!"
    exit 1
fi

# Set default values to parameters that haven't been specified
if [ "$sleep_interval" == "" ]; then sleep_interval="1m"; fi
if [ "$start_step" == "" ]; then start_step=00; fi
if [ "$label" == "" ]; then label=$run_dir_name; fi

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

# Make sure that force-step is greater than or equal to start-step
if [ "$force_step" != "" ] && [ "$force_step" -lt "$start_step" ]
then
  if [ "$start_step" -gt 00 ]
  then
    echo "Error: Value of --force-step $force_step must not be less than that of --start-step $start_step. Exiting..."
    exit 1
  else
    echo "Error: Value of --force-step $force_step must not be less than $start_step. Exiting..."
    exit 1
  fi
fi

# Default in pipeline_type if this is a rerun
if [ "$pipeline_type" == "" -a -e "`ls $run_dir/progress-*.txt 2>/dev/null`" ]
then pipeline_type=`basename $run_dir/progress-*.txt .txt | sed 's/^progress-//'`
fi

# Read pipeline config and 
# Set file where pipeline progress will be written to
if [ "$pipeline_type" != "" ]
then
  which_config.sh app $pipeline_type
  exit_status=$?
  if [ "$exit_status" -eq 0 ]; then . import_config.sh app $pipeline_type
  else echo "Note: Ignoring absent $pipeline_type config."
  fi
else
  pipeline_type=pipeline
fi
pipeline_progress_file=progress-${pipeline_type}.txt

# Change to the run directory
cd $run_dir
exit_status=$?
if [ "$exit_status" -ne 0 ]; then exit $exit_status; fi

if [ "$AFC_STEPS" == "" ] 
then 
  AFC_STEPS=$(ls $run_dir/*.sh | sed -re "s#$run_dir/##" | grep -e "^[0-9]\+.\.sh" | sed -re "s/.\.sh//" | uniq)
fi 

preprocess

# Run each step of the pipeline
echo "Starting $pipeline_type pipeline..."
for step in $AFC_STEPS
do
  # If current step is less than the step from which pipeline should be 
  # started, skip to next iteration
  if [ "$step" -lt "$start_step" ]
  then 
    echo "Skipping step $step because it's less than the start-step $start_step"
    continue 
  fi
  # If current step is greater than the max step the pipeline should run,
  # finish.
  if [ "$max_step" != "" ] && [ "$step" -gt "$max_step" ]
  then
    notify_success "Finished max step $max_step"
    exit 
  fi

  # Run step
  if [ "$step" == "$force_step" ]
  then 
    force='force'
  else
    force= 
  fi
  run_step.sh $step $force >$tmp_file 2>&1
  exit_status=$?
  cat $tmp_file | log -
  if [ "$exit_status" -ne 0 ]
  then failure "run_step.sh exited with code $exit_status; perhaps there was a QC failure" $exit_status
  fi
  
  # Find the LSF job ID of the job that was submitted
  lsf_job_id=`get_lsf_job_id.sh $tmp_file`
  exit_status=$?
  if [ "$exit_status" -ne 0 ]
  then
    failure "get_lsf_job_id.sh exited with code $exit_status; either get_lsf_job_id.sh failed or job was not submitted to cluster."
  elif [ "$lsf_job_id" == "" ]
  then
    log "Status: Started (LSF job ID: n/a)"
    log "Status: Completed"
    continue
  fi
  lsf_job_id=`echo $lsf_job_id | tr "\n" " " | sed 's/\s$//'`
  log "Status: Started (LSF job ID: $lsf_job_id)"
  
  # Check whether job has completed, and wait till it completes
  failed_bjobs_out="$run_dir/failed_bjobs_out_"`echo "$lsf_job_id" | tr ' ' _`".txt"
  check_lsf_job_status.sh "$lsf_job_id" $sleep_interval $failed_bjobs_out
  exit_status=$?
  if [ "$exit_status" -eq 2 ]
  then
    log "Error: At least one job failed (based on job status):"
    cat $failed_bjobs_out | log -
    failure "Job scheduler reported job failure" $exit_status
  else
    # Go through the stdout files and make sure none of the jobs failed. This 
    # is necessary because sometimes it's possible to miss detecting a failed 
    # job by just using bjobs.
    # Before checking the log files, make sure they are all there. This is 
    # necessary because sometimes the file isn't fully written right away.
    log_dir=`echo $run_dir/file*-$step*/logs`
    exp_num_jobs=`cat $run_dir/cmds-$step.sh | wc -l`
    if [ "$exp_num_jobs" == 0 ]
    then
      log "Could not determine job count; sleeping 30s before checking logs"
      sleep 30
    else
      # Maximum number of sleeps
      remaining_sleeps=10
      curr_sleep=2
      # Wait for all files to exist, and get list of files
      log_list_file=`mktemp`
      while true
      do
        # Build list of files
        for i in $lsf_job_id
        do ( cd $log_dir; echo *out.${i}.*.txt ) >> $log_list_file
        done
        # If it's too short, then sleep and iterate
        if [ `cat $log_list_file | wc -w` -lt $exp_num_jobs ]
        then
          log "Expecting $exp_num_jobs, but found only `cat $log_list_file | wc -w` log files"
          let --remaining_sleeps
          let 'curr_sleep *= 2'
          if [ "$remaining_sleeps" -lt 0 ]; then break; fi
          log "Will check again in $curr_sleep"
          cat /dev/null > $log_list_file
          sleep $curr_sleep
        # Otherwise break out
        else
          log "All $exp_num_jobs logs exist"
          break
        fi
      done
      # Reset sleeps, unless we had already run out
      if [ "$remaining_sleeps" -ge 0 ]
      then
        remaining_sleeps=10
        curr_sleep=2
      fi
      # Now, wait for all of the files to be non-empty
      still_empty_log_file=`mktemp`
      while [ "$remaining_sleeps" -ge 0 ]
      do
        # Get list of empty files
        for file in `cat $log_list_file`
        do if [ ! -s "$log_dir/$file" ]; then echo $file; fi
        done > $still_empty_log_file
        # If there are some, then sleep and iterate
        if [ -s $still_empty_log_file ]
        then
          log "There are still `cat $still_empty_log_file | wc -w` empty log files"
          cat $still_empty_log_file | log -
          let --remaining_sleeps
          let 'curr_sleep *= 2'
          if [ "$remaining_sleeps" -lt 0 ]; then break; fi
          log "Will check again in $curr_sleep"
          cp $still_empty_log_file $log_list_file
          sleep $curr_sleep
        # Otherwise break out
        else
          log "All logs are non-empty"
          break
        fi
      done
      # Now, either all files exist, or we gave up trying
      if [ "$remaining_sleeps" -lt 0 ]
      then failure "Timed out waiting for log files" 3
      fi
    fi
    
    echo "Checking stdout in $log_dir to see if any jobs failed..."
    for i in $lsf_job_id;
    do
      result="`grep -L "Success" $log_dir/*out.${i}.*.txt`"
      if [ ! -z "$result" ]
      then
        log "Error: At least one job failed (based on log):"
        echo "$result" | log -
        failure "Job failure detected in logs" 2
      fi
    done
    echo "...No jobs failed in this step!"
  fi 
  log "Status: Completed"
done

# Finally, indicate in the pipeline-progress-file that the pipeline has completed
log "Status: Pipeline $pipeline_type completed!"
completion
