#!/bin/bash
# Writes a script to run swap_symlink.sh on a set of symlink paths
#
# If you specify one or more directories as arguments, then commands will be
# written for each symlink directly contained within those directories
#
# Parameters:
# $1,... = paths to symlinks or directories

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
PATHS="$@"

# Iterate over paths and write commands
for path in $PATHS
do
  # Get symlink list for the path
  if [ -d $path ]
  then
    symlinks=
    for file in $path/*
    do if [ -h $file ]; then symlinks="$symlinks $file"; fi
    done
  else
    symlinks=$path
  fi
  
  # Write comands
  if [ "$symlinks" == "" ]; then continue; fi
  for symlink in $symlinks
  do echo swap_symlink.sh $symlink
  done
done
