#!/bin/bash
# Generic wrapper for the varying array submission types. 
# Imports the correct one to use from the config util-scripts

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

. import_config.sh app util-scripts AFC_SCRIPT

$AFC_SCRIPT "$@" 
