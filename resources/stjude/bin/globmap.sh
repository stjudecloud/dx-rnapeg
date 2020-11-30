#!/bin/bash
# BASH implementation of Perl's globmap
#
# This first parameter is a glob, and the second is an output spec.  The return
# is tab-delimited text where the first column is a file path matching the glob,
# and the second is an output path built using the input and spec.
#
# For more information, see perldoc for File::GlobMapper
#
# For example, if foo contains files a.tar.gz, b.txt.gz, c.zip, then:
#
# globmap.sh foo/*.gz bar/#1
#
# returns
# 
# foo/a.tar.gz  bar/a.tar
# foo/b.txt.gz  bar/b.txt
#
# Parameters:
# $1 = glob pattern (IMPORTANT: be sure to quote it to avoid shell expansion)
# $2 = output spec, using #n as capture groups

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Error out if 1 or 3+ parameters were sent
if [ "$#" != 2 ]
then
  about.sh $0
  echo
  echo "Must use exactly two parameters." >&2
  echo "Be sure to quote the glob!" >&2
  exit 1
fi

# Get parameters
GLOB="$1"
SPEC="$2"

regex=`echo "$GLOB" | sed -e 's#\?#([^/])#g' -e 's#\*#([^/]*)#g'`
replace=`echo "$SPEC" | sed -r 's/#([0-9])/\\\1/g'`
sedcmd="s/^$regex\$/\\0\\t$replace/"

echo "regex=$regex"
echo "replace=$replace"
echo "sedcmd=$sedcmd"

ls -d $GLOB | sed -r "$sedcmd"