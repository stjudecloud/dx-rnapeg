#!/bin/bash
# Parses a sample pair string into its sample names
#
# The two names are written to stdout with a space between
#
# Usage (assumes sample pair string to parse is in $samplepair):
# read case control < <(sn_string_to_pair.sh $samplepair)
#
# There is no real error checking here
#
# $1 = sample pair string

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
PAIR=$1

# Parse
echo $PAIR | gawk -F _ '
{
  # Assemble case sample and control string
  caseSample = $1 "_" $2
  
  # Assemble control sample
  if(length($3) <= 2) {
    controlSample = $1 "_" $3
  }
  else {
    controlSample = $3 "_" $4
  }
  
  # Return samples
  print caseSample, controlSample
}'
