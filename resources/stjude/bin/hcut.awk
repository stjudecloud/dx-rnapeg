#!/bin/gawk -f
# Re-implementation of cut that allows cutting by header name.
#
# The field separator you specify in the parameters will override what you tell
# awk on the command line, and the default is tab.  There is no option for an
# all-whitespace separator.
#
# The input must have a single header line.
#
# The output has NO HEADER by default.
#
# Parameters:
# fieldsep: field separator (default tab)
# fields: comma-separated list of field names or numbers to emit
# header: 0, blank, or omitted to write no header, otherwise write a header
# lenient: true value to allow invalid field names in fields
# mode: require, optional, remove (default require)
BEGIN {
  # Show usage information if some parameters were missing
  if(fields == "") {
    print "Missing required parameters, will attempt to print usage."
    system("about.sh `which hcut.awk`")
	exit 1
  }
  
  # Initialize
  if(!fieldsep) fieldsep = "\t"
  FS = fieldsep
  OFS = FS
  if(!mode) mode = "require"
  if (mode == "optional") lenient=1
  split("", header2col, "")
}
function convertToColNumber(nameOrNumber) {
  if(0 + nameOrNumber < 1) {
    col = header2col[nameOrNumber]
    if(col == 0 && !lenient) {
      print "Invalid field name: " nameOrNumber > "/dev/stderr"
      exit 1
	}
    return col
  }
  else {
    return nameOrNumber
  }
}
function convertToArrayOfColNumbers(list, array) {
  lengthof_arr = split(list, arr, ",")  
  for(i = 1; i <= lengthof_arr; i++) array[i] = convertToColNumber(arr[i])
  return lengthof_arr
}
function printCutLine() {
  j = 1
  for(i = 1; i <= numFields; i++) {
    fieldNum = fieldNums[i]
    if(fieldNum == 0) continue
    if(j > 1) printf fieldsep
    printf $fieldNum
    j++
  } 
  print ""
}
NR == 1 {
  # Read header
  numcols = NF
  for(i = 1; i <= numcols; i++){
    header2col[$i] = i
  }
 
  # Get column numbers to print

  if (mode != "remove"){
    numFields = convertToArrayOfColNumbers(fields, fieldNums)
  }
  else{
    numRemoveFields = convertToArrayOfColNumbers(fields, fieldRemoveNums)
	for(i = 1; i <= numRemoveFields; i++) removeSet[fieldRemoveNums[i]] = 1;
	numFields = 0
	for(i = 1; i <= numcols; i++) {
		if(!removeSet[i]) {
			fieldNums[++numFields] = i
		}
	}
  } 
  # Print it, if appropriate
  if(header) printCutLine()
  next
}
{
  printCutLine()
}
