#!/bin/bash
# Renames a file to put the file's modified time in its filename just before
# the extension
#
# By default this will place the suffix before the last dot.  If your file has
# multiple extensions, then you can specify the number of true extension parts
# in the second parameter:
# 
# timestamp.sh foo.bar.txt     -> foo.bar_TIMESTAMP.txt
# timestamp.sh foo.bar.txt 2   -> foo_TIMESTAMP.bar.txt
#
# If you specify 0 extension parts, then it will be appended to the end of the
# filename.
#
# If you specify more extension parts than there actually are, then it will be
# inserted before the first dot, e.g.
#
# timestamp.sh foo.bar.txt 9   -> foo_TIMESTAMP.bar.txt
#
# This will NOT rename a file over an existing file (except under a race
# condition)
#
# If the move is successful, then the new filename will be written to stdout,
# and a 0 exitcode is returned.  If the move fails for any reason, then
# nothing is written to stdout, and a non-0 exit code is returned.
#
# $1 = file to rename
# $2 = (optional) how many true extension parts there are (default no limit)

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
FILE=$1
PARTS=$2
if [ ! $PARTS ]; then PARTS=1; fi

# Make sure file exists
if [ ! -e $FILE ]
then echo "File not found: $FILE" >&2; exit 1
fi

# Get and format the timestamp
ts=`stat -c '%y' $FILE`
ts=`date -d "$ts" +%Y%m%d%H%M`

# Get old path, name, and new name
DIR=`dirname $FILE`
oldname=`basename $FILE`
newname=`echo $oldname | awk -v parts=$PARTS -v ts=$ts '
BEGIN { FS="." }
{
  last = NF - parts
  if(last < 1) last = 1
  newname = ""
  for(i = 1; i <= last; i++) {
    if(newname) newname = newname "."
    newname = newname $i
  }
  newname = newname "_" ts
  for(i = last + 1; i <= NF; i++) {
    newname = newname "." $i
  }
  print newname
}'`

newfile=$DIR/$newname
if [ -e $newfile ]
then echo "Will not overwrite existing file $newfile" >&2; exit 1
fi

if mv -nv $FILE $newfile >&2
then
  if [ -e $FILE ]
  then echo "Move failed" >&2; exit 1
  else echo $newfile
  fi
else exit $?
fi
