#!/bin/bash
# Runs the command found on a single line of an input file.  Use this to run
# each line of a script as a separate task under a single SGE array job.
#
# Usage:
# qsub -t 1-`cat script.sh | wc -l` sge_array_for_cmdfile.sh script.sh
#
# $1 = file with one task on each line

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

if [ "$SGE_TASK_ID" == "" ]
then echo "No SGE_TASK_ID; must do qsub -t to use this script" >&2; exit 1
fi
cmd=`head -n $SGE_TASK_ID $1 | tail -n 1`
eval $cmd
