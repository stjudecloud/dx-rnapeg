#!/bin/sh
# Prints a file's header and absolute real path (w/ special handling for Java
# class files).
# The header is the whole part up to the first non-blank line that is not a
# comment.
#
# If the file is a java class file, then it runs it via java_sj_compbio.sh with
# no arguments; by convention, it should print usage instructions.
#
# $1 = the file

file=$1

# Check for java class file
if [ "${file: -6}" == ".class" ]
then
  file=`readlink -f $file`
  classname=`basename $file | sed 's/.class//'`
  path=`dirname $file`
  while [ -n "$path" ]
  do
    part=`basename $path`
    classname="$part.$classname"
    if [ "$part" == "org" ]; then break; fi
    path=`dirname $path`
  done
  java.sh $classname
  echo
  echo java.sh $classname
else
  if [ ! -f "$file" ]; then file=`which $file`; fi
  readlink -f $file
  echo
  sed -r -e '/^($|[^#])/,$d' -e '/^#!/d' -e 's/^# ?//' $file
fi
