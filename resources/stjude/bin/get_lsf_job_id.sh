#!/bin/bash
# Extracts LSF job ID from the given file containing bsub output and prints it
# to stdout.
#
# $1 = File containing bsub output from which LSF job ID is to be extracted

# Show usage information if no parameters were sent
if [ $# -lt 1 ]; then about.sh $0; exit 1; fi

input_file=$1
cat $input_file | grep -e "is submitted to .*queue" | cut -d " " -f2 | tr -d "<" | tr -d ">"
