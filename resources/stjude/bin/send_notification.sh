#!/bin/bash
# Sends a notification to email addresses read from the notifications config file.
#
# $1 = notification name
# $2 = sequencing type
# $3 = project
# $4 = subject
# $5 = file containing message body; omit to read from stdin

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
NAME=$1
TYPE=$2
PROJECT=$3
SUBJECT="$4"
BODY_FILE=$5

# Default to stdin if no body file
if [ "$BODY_FILE" == "" ]
then BODY_FILE=/dev/stdin
fi

# Get email addresses
addrs=`mktemp`
awk -v name=$NAME -v type=$TYPE -v project=$PROJECT '
  /^#/ { next }
  $1 == name
      && ($2 == "*" || $2 == type)
      && ($3 == "*" || $3 == project) {
    print $4
  }' | sort | uniq > $addrs

# Send email
if [ -s "$addrs" ]
then
  cat $BODY_FILE | Mail -s "$SUBJECT" `cat $addrs`
  if [ "$?" == 0 ]
  then echo "Mail sent successfully" >&2
  else echo "Failure in Mail command" >&2; exit 1
  fi
else
  echo "Nobody is subscribed to this notification" >&2
fi