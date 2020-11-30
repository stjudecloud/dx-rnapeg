#!/bin/bash
# Builds a sample name or barcode from its parts
#
# Sample name/barcode is written to stdout.
#
# Usage (assumes inputs are in variables):
# samplename=`sn_build.sh $subject $disease $type $index $sid`
#
# On build error, nothing is written to stdout, and messages are written to
# stderr
#
# $1 = subject
# $2 = disease
# $3 = type (single letter code)
# $4 = index (optional; default is 1; required if secondary ID is specified)
# $5 = secondary ID (optional; iff included, output is a barcode)

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters (awk will parse)
INPUT="$@"

result=`echo "$INPUT" | gawk --re-interval '
# Init
BEGIN {
  sample = ""
  numtypes = split("A C D G M O R X", types);
  for(ord = 65; ord < 91; ord++) {
    ords[sprintf("%c", ord)] = ord;
  }
  error = ""
}
# Parse input, and make it available to END block
{
  subject = $1
  disease = $2
  type = $3
  idx = $4
  sid = $5

  # Validate subject
  if(! (subject ~ /^SJ/)) {
    error = "Invalid subject " subject
    exit
  }
  # Validate point
  typeok = 0
  for(i in types) if(type == types[i]) { typeok = 1; break }
  if(!typeok) {
    error = "Invalid type " type
    exit
  }
  
  ## check if index is larger than zero 
  if(idx <= 0) {
     error = "Index should be larger than zero"
   }

  # Default index
  if(idx == "") idx = 1;
  
  # Suffix
  if(sid == "") suffix = ""
  else suffix = "-" sid
}
# V1
$1 ~ /^SJ[A-Za-z0-9]*[A-Za-z][0-9]{3,4}$/ {
  # Check subject-disease consistency
  if(! (subject ~ ("SJ" disease "[0-9]{3}") ) ) {
    error = "Inconsistent subject and disease: " subject ", " disease
    exit
  }
  
  # Process point
  pointord = ords[type] + idx - 1;
  pointcode = sprintf("%c", pointord)
  
  # Build sample
  sample = subject "_" pointcode suffix
  
  # Done, go to end
  exit
}
# Default: V2
{
  sample = "SJ" disease substr(subject, 3) "_" type idx suffix
}
# Output
END {
  # Handle error
  if(error) {
    print "ERROR " error
    exit
  }
  
  # Output
  print sample
}
'`

if [ "${result:0:5}" == "ERROR" ]
then echo $result >&2
else echo $result
fi
