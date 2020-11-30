#!/bin/bash
# Builds a glob pattern to match files for a subject.
#
# Glob pattern is written to stdout
#
# Usage (assumes inputs are in variables):
# glob=`sn_subj_glob.sh $subject`
#
# On error, nothing is written to stdout, and messages are written to stderr
#
# Note: the glob pattern always begins with SJ and ends with *
#
# $1 = subject

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters (awk will parse)
INPUT="$@"

result=`echo "$INPUT" | gawk --re-interval '
# Init
BEGIN {
  glob = ""
  error = ""
}
# Parse input, and make it available to END block
{
  subject = $1
  #disease = $2

  # Validate subject
  if(! (subject ~ /^SJ/)) {
    error = "Invalid subject " subject
    exit
  }
}
# V1
$1 ~ /^SJ[A-Za-z]+[0-9]{3,4}$/ {
  # Check subject-disease consistency
  #if(! (subject ~ ("SJ" disease "[0-9]{3}") ) ) {
  #  error = "Inconsistent subject and disease: " subject ", " disease
  #  exit
  #}
  
  # Build glob
  glob = subject "_*"
  
  # Done, go to end
  exit
}
# Default: V2
{
  glob = "SJ*" substr(subject, 3) "_*"
}
# Output
END {
  # Handle error
  if(error) {
    print "ERROR " error
    exit
  }
  
  # Output
  print glob
}
'`

if [ "${result:0:5}" == "ERROR" ]
then echo $result >&2
else echo "$result"
fi
