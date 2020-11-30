#!/bin/bash
# Runs multiple pipelines in the given order
#
# Parameters:
# -p STR  List of pipelines. Each pipeline should be of the form 
#          "<pipeline type>:<pipeline run dir>" (double-quotes are required). 
#          Multiple pipelines can be specified. The pipelines that are to be run
#          together in parallel should be separated by a semi-colon, and the 
#          pipelines that are to be run in sequence should be separated by a space. 
#          For example: -p "coverage:/path/to/run/dir;snp:/path/to/run/dir coverage-post:/path/to/run/dir;crest:/path/to/run/dir cnv:/path/to/run/dir"
#          In this example, coverage and snp pipelines will run in parallel. 
#          Next, once the coverage and snp pipelines are complete, coverage-post 
#          and crest pipelines will run in parallel. Once those are finished, the 
#          cnv pipeline will run.
# -n STR  (Optional) Comma-separated email addresses that should be notified when 
#           all the pipelines complete or if any pipeline fails 
#           (Ex: john.doe@stjude.org,jane.smith@stjude.org).
#           The value of the environment variable USER_EMAIL is appended to this 
#           list. 
# -N      (Optional) Individual pipeline notifications will be sent for 
#           job failures only and not for successful job completions.
# -h      (Optional) Show usage.

# Show usage information if fewer parameters were sent
if [ "$#" -lt 1 ]; then about.sh $0; exit 1; fi

# Notify people of pipeline completion or failure
# $1 = Recipient email addresses, comma-separated 
# $2 = Subject
# $3 = From
# $4 = Message
send_email() {
  local from=$1
  local to=$2
  local subject=$3
  local message=$4
  echo "$message" | mail -r "$from" -s "$subject" "$to"
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
pipeline_list=
notify_recipients=
failure_notification_only=
while getopts "hp:n:N" OPTION
do
  case $OPTION in
    h)
      about.sh $0
      exit 1
      ;;
    p)
      pipeline_list=$OPTARG
      ;;
    n)
      notify_recipients=$OPTARG
      ;;
    N)
      failure_notification_only="failure_notification_only"
      ;;
  esac
done

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

# Set notification param string
notification_param_str=
if [ ! "$failure_notification_only" ]
then
  notification_param_str="-n '$notify_recipients'"
else
  notification_param_str="-n '$notify_recipients' -N"
fi

# Run pipelines in proper order
for pipeline_set in $pipeline_list
do
  pipelines=`echo $pipeline_set | sed 's/;/ /g'`
  process_id_list=
  for pipeline in $pipelines
  do
    read part1 part2 < <(echo $pipeline | tr : ' ')
    if [ "$part2" == "" ]
    then
      run_dir="$part1"
      pipeline_type=`basename $run_dir`
      cmd="run_pipeline.sh -d $run_dir $notification_param_str"
    else
      run_dir="$part2"
      pipeline_type="$part1"
      cmd="run_pipeline.sh -d $run_dir -t $pipeline_type $notification_param_str"
    fi
    echo "Submitting $pipeline_type pipeline having run dir ${run_dir}..."
    echo "$cmd"
    eval "$cmd" &
    process_id_list="$process_id_list $!"
  done
  pipeline_failed=
  for i in $process_id_list
  do
    wait $i
    exit_code=$?
    if [ "$exit_code" -ne 0 ]; then pipeline_failed="yes"; fi
  done
  if [ "$pipeline_failed" == "yes" ]
  then
    echo "Error: At least one pipeline failed!"
    subject="Failure in one of the pipelines"
    message="At least one pipeline in the pipeline list (${pipelines}) failed"
    notify "$subject" "$message"
    exit 1
  fi
  echo "Completed pipelines $pipelines"
done

# Print pipeline-completed status
echo "Your pipelines have successfully completed!"
subject="Your pipelines have completed!"
message="Your pipelines $pipeline_list have completed!"
notify "$subject" "$message"
