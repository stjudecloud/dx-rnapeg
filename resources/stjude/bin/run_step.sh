#!/bin/bash
# Takes a step number and runs all of the substeps associated with that step
#
# Parameters: 
# $1 = step number (01, 02, etc.)

# Show usage information if no parameters were sent
if [ "$#" -lt 1 ]; then about.sh $0; exit 1; fi

. steplib.sh
CWD=$( pwd ) 
set_step_script_dir $CWD
STEP=$1
FORCE=$2

run_step $STEP $FORCE 
