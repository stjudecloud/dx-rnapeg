#!/bin/bash
# Finds the config file for a configuration name and echoes it.
#
# This searches for a file named $1.config.txt in the category configuration
# directory under the root configuration directory.
# The root configuration directory used will be:
# The directory in the SJ_CONFIGS environment variable
#
# Usage: configfile=`which_config.sh mycategory myconfig`
#    OR: configrootdir=`which_config.sh .`
#    OR: configdir=`which_config.sh mycategory`
#
# $1 = category name, or . for the root category
# $2 = config name or abs path to file (starts with /), or omitted to just
#      report the config directory

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
CATEGORY=$1
CONFIG=$2

# Check for absolute path case
if [ "${CONFIG:0:1}" == "/" ]
then
  if [ -f "$CONFIG" ]
  then
    readlink -f $CONFIG
    exit
  else
    echo "No config file at given absolute path: $CONFIG" >&2
    exit 1
  fi
fi

# Find root config dir
# Use SJ_CONFIGS environment variable
if [ $SJ_CONFIGS ]
then dir=$SJ_CONFIGS
else
  about.sh $0
  echo -e "CONFIGURATION DIRECTORY COULD NOT BE FOUND\nMAKE SURE ENV VARIABLE SJ_CONFIGS IS SET" >&2
  exit 1
fi

# Go down into the category dir
dir=$dir/$CATEGORY

# Report directory if config was omitted
if [ ! $CONFIG ]; then echo $dir; exit; fi

# Find config file
config_file=$dir/$CONFIG.config.txt

# If it was found, then write it out, and we are done
if [ -f "$config_file" ]
then
  readlink -f $config_file
  exit
fi

# Otherwise, look for alias file, and see if config is an alias
alias_file=$dir/aliases.txt
if [ -f $alias_file ]
then newconfig=`awk -v alias=$CONFIG '$1 == alias { print $2; exit }' $alias_file`
else newconfig=
fi

# If an alias entry was found, then try again under the real name
if [ $newconfig ]
then
  echo "Config loaded using alias $CONFIG for $newconfig" >&2
  config_file=$dir/$newconfig.config.txt

  # Again, if it was found, then write it out, and we are done
  if [ -f "$config_file" ]
  then
    readlink -f $config_file
    exit
  fi
fi

# Finally, throw error
echo "Could not find config for $CONFIG at $config_file" >&2
exit 64
