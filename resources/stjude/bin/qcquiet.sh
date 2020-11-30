#!/bin/bash
# Runs a QC script and saves stdout and stderr iff the result is not a
# pass; writes last-line result to stdout for reading by caller, and exits
# with a non-zero exit code if result is not PASS.
#
# Uses hard-coded temp files .lastqc.err.txt and .lastqc.out.txt in the
# current working dir.  These files should be deleted on exit.
#
# If the script passes, then no files are created.
#
# See qclib.sh
#
# $1 = output prefix (.out.txt and .err.txt will be appended)
# $2,... = qc command and params

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
PREFIX=$1
shift
CMD="$@"

result=`$CMD 2> .lastqc.err.txt | tee .lastqc.out.txt | tail -n 1`
echo "$result [`basename $PREFIX`]"
if [ "${result:0:4}" == "PASS" ]
then
  rm .lastqc.out.txt .lastqc.err.txt
  exit 0
else
  mv .lastqc.out.txt $PREFIX.out.txt
  mv .lastqc.err.txt $PREFIX.err.txt
  exit 1
fi
