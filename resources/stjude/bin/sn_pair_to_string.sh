#!/bin/bash
# Creates a sample pair string from a pair of sample names
#
# The sample pair string is written to stdout
#
# Usage (assumes sample names are in $case and $control):
# pair=`sn_pair_to_string.sh $case $control`
#
# There is no real error checking here
#
# $1 = sample pair string

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters (awk will parse)
INPUT="$@"

# Parse (strip off SIDs before sending to awk)
echo "$INPUT" | sed -r 's/-[^ ]*( |$)/ /g' | gawk '
{
  # Check for match up to first underscore, which indicates short form may be
  # used
  pos = index($1, "_")
  if(substr($1, 1, pos) == substr($2, 1, pos)) {
    print $1 "_" substr($2, pos + 1)
  }
  else {
    print $1 "_" $2
  }
}'
