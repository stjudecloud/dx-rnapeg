#!/bin/bash
#
# DEPRECATED July 6, 2016. Use java.sh instead as it encorporates this
# functionality by default.
#
# $1 = directory to use as java.io.tmpdir (if "-", then it calls mktemp -d)
# $2... = arguments sent to java.sh (and then to java itself)

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get the temp dir to use
#JAVA_TMPDIR=$1
shift

# If TMPDIR == "-", then default to mktemp -d
#if [ "$JAVA_TMPDIR" == "-" ]
#then JAVA_TMPDIR=`mktemp -d`
#fi

# Delegate to the normal java.sh, with the tmpdir param added
echo
echo "################################################################"
echo "# java-settmp.sh is now deprecated. Please use java.sh instead #"
echo "################################################################"
echo 
java.sh "$@"
