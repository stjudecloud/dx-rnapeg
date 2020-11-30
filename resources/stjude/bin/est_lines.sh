#!/bin/sh
# Estimates the number of lines in a file by sampling the head of the file
#
# Also works for gzipped files
#
# $1 = file whose size to estimate
# $2 = (optional) number of bytes to sample (default 900000)

# Get arguments
FILE=$1
SAMPLE_SIZE=$2

if [ "$SAMPLE_SIZE" == "" ]; then SAMPLE_SIZE=900000; fi

FILE=`readlink -f $FILE`
if [ ! -e "$FILE" ]
then 
  exit 1
fi
filesize=`stat -c '%s' $FILE`

if [ "$filesize" -le $SAMPLE_SIZE ]
then samplecmd="cat"; SAMPLE_SIZE=$filesize
else samplecmd="head -c $SAMPLE_SIZE"
fi

if file `readlink -f $FILE` | grep -q gzip
then gunzipcmd="gunzip -c"
else gunzipcmd="cat"
fi

#echo "filesize=$filesize SAMPLE_SIZE=$SAMPLE_SIZE" >&2
eval $samplecmd $FILE | eval $gunzipcmd 2>/dev/null | wc -l | \
  ( read lines ; let "est = $filesize * $lines / $SAMPLE_SIZE"; echo $est )
