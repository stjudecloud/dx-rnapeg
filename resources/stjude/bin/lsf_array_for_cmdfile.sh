#!/bin/bash
# Runs the command found on a single line of an input file.  Use this to run
# each line of a script as a separate task under a single LSF array job.
#
# Usage:
# bsub -J NAME[1-`cat script.sh | wc -l`] lsf_array_for_cmdfile.sh script.sh
#
# $1 = file with one task on each line

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

if [ "$LSB_JOBINDEX" == "" ]
then echo "No LSB_JOBINDEX; must do bsub as job array to use this script" >&2; exit 1
fi
cmd=`head -n $LSB_JOBINDEX $1 | tail -n 1`
eval "$cmd"
