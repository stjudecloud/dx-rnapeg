#!/bin/bash
# Parses a sample name into its parts
#
# Parts are written in a space-delimited list to stdout.
#
# Usage (assumes sample name to parse is in $samplename):
# read subject disease type index sid < <(sn_parse.sh $samplename)
#
# On parse error, nothing is written to stdout, and messages are written to
# stderr
#
# $1 = sample name

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
BARCODE=$1

result=`echo $BARCODE | gawk --re-interval '
# Init
BEGIN {
  input = ""
  subject = ""
  disease = ""
  type = ""
  idx = 0
  sid = ""
  numtypes = split("A C D G M O R X", types);
  for(ord = 65; ord < 91; ord++) {
    ords[sprintf("%c", ord)] = ord;
  }
  error = ""
}
# Make input available to END block
{ input = $0 }
# V1
/^SJ[A-Za-z0-9]*[A-Za-z][0-9]{3,4}_/ {
  # Main regex match to capture groups
  if(!match($0, /^SJ([A-Za-z0-9]+)([0-9]{3})_([A-Z])(-(.+))?$/, parts)) {
    error = "Could not parse " $0 " using format SJdddnnn_t"
    exit
  }
  
  # Get parts
  disease = parts[1];
  num = parts[2];
  pointcode = parts[3];
  sid = parts[5];
  
  # Process point
  pointord = ords[pointcode]
  lastord = 0
  for(i = 1; i <= numtypes; i++) {
    code = types[i]
    ord = ords[code]
    if(ord > pointord) break
    lastord = ord
  }
  type = sprintf("%c", lastord)
  idx = pointord - lastord + 1
  
  # Build subject
  subject = "SJ" disease num
  
  # Done, go to end
  exit
}
# Default: V2
{
  # Main regex match to capture groups
  if(!match($0, /^SJ([A-Za-z0-9]+)([0-9]{6})_([A-Z][0-9]+)(-(.+))?$/, parts)) {
    error = "Could not parse " $0 " using format SJdddnnnnnn_ti"
    exit
  }
  
  # Get parts
  disease = parts[1];
  num = parts[2];
  pointcode = parts[3];
  sid = parts[5];

  # Process point  
  type = substr(pointcode, 1, 1)

  ## check if the type entered is valid 
  validType = 0
  for(i = 1; i <= numtypes; i++) {
     if (type == types[i]) {
         validType = 1
     }
  }
  if (!validType) {
     error = "Invalid type, valid types for v2 are A, C, D, G, M, O, R and X"
  }
  
  idx = substr(pointcode, 2)
  
  if (idx == 0){
     if (error) {
        error = error "; Invalid index, index should be non-zero"
     } 

     else {
        error = "Invalid index, index should be non-zero"
     }
  }

  # Build subject
  subject = "SJ" num;
}
# Output
END {
  # Handle error
  if(error) {
    print "ERROR " error
    exit
  }
  
  # Validate disease
  if(disease ~ /[0-9][0-9]$/) {
    print "ERROR Could not parse " input ", got bad disease code " disease
    exit
  }
  
  # Output
  print subject, disease, type, idx, sid
}
'`

if [ "${result:0:5}" == "ERROR" ]
then echo $result >&2
else echo $result
fi
