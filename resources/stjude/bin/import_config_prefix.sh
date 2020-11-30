#!/bin/bash
# Imports one or more values from a config file based on prefix; this must be 
# sourced to have any effect.
#
# Example:
# . import_config_prefix.sh genome GRCh37-lite REFSEQ_NM_ -t
#
# Note that this will exit on error; since you are sourcing it that means that
# it has the ability to abort the calling script.
#
# $1 = category name
# $2 = config name 
# $3 = prefix to import
# $4 = (optional) -t to trim the prefix off the variable names

# Find the config file
_CONFIG_FILE=`which_config.sh $1 $2`
if [ ! -f $_CONFIG_FILE ]; then echo "No config found; exiting"; exit 65; fi

# Get other parameters
_PREFIX=$3
_TRIM=$4

# Load based on prefix
_search="^$_PREFIX[^\t]*\t"
while read _assign
do
  if [ "$_TRIM" == "-t" ]; then _assign=`echo "$_assign" | sed "s/^$_PREFIX//"`; fi
  a=$( echo  "$_assign" | cut -f 1 ) 
  b=$( echo  "$_assign" | cut -f 2 )
  read ${a} <<IN
$b
IN
done <<< "`grep -P "$_search" $_CONFIG_FILE`"
