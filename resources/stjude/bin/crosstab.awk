#!/bin/gawk -f
# Takes tabular text input data and writes out cross-tabulated data, similar to
# a PivotTable in Excel.
#
# Input must be sorted by the rowcols, and must be a file because multiple
# passes are made.
#
# The field separator you specify in the parameters will override what you tell
# awk on the command line, and the default is tab.  There is no option for an
# all-whitespace separator.
#
# The output will always have a header, so it is suggested you specify one on
# input.
#
# Parameters:
# fieldsep: field separator (default tab)
# header: 0, blank, or omitted if no header, otherwise there is a header
# numcols: number of columns (req'd and used iff no header)
# rowcols: row-defining columns (comma-separated list of names or numbers)
# colcols: col-defining columns (comma-separated list of names or numbers)
# xtabcol: crosstab column name or number (1-based)
# xtabfun: function (min, max, sum, first, last, or catDELIM), default last
# morefun: function for extra columns (not row, col, or xtab)
# morefuns: list of colnum:function for "more" column function overrides
BEGIN {
  # Show usage information if some parameters were missing
  if(xtabcol == "" || rowcols == "" || colcols == "" || (!header && !numcols)) {
    print "Missing required parameters, will attempt to print usage."
    system("about.sh `which crosstab.awk`")
	  exit 1
  }
  
  # Process morefuns
  split("", morefunsarr, "")
  if(morefuns) {
    split(morefuns, parts, ",")
	for(i = 1; i <= length(parts); i++) {
	  split(parts[i], pair, "=")
	  morefunsarr[pair[1]] = pair[2]
	}
  }
  
  # Initialize
  if(!fieldsep) fieldsep = "\t"
  FS = fieldsep
  OFS = FS
  split("", header2col, "")
}
function convertToColNumber(nameOrNumber) {
  if(0 + nameOrNumber < 1) return header2col[nameOrNumber]
  else return nameOrNumber
}
function convertToArrayOfColNumbers(list, array) {
  lengthof_arr = split(list, arr, ",")
  for(i = 1; i <= lengthof_arr; i++) array[i] = convertToColNumber(arr[i])
}
function getKey(row, cols) {
  key = row[cols[1]]
  for(i = 2; i <= length(cols); i++) key = key "," row[cols[i]]
  return key
}
NR == 1 {
  # Read header if there is one
  if(header) {
    numcols = NF
    for(i = 1; i <= numcols; i++) header2col[$i] = i
	header = $0
  }

  # Convert *cols to arrays of numbers
  convertToArrayOfColNumbers(rowcols, rowcolsarr)
  convertToArrayOfColNumbers(colcols, colcolsarr)
  xtabcol = convertToColNumber(xtabcol)
  
  # Remember which columns are colcols, as they should be suppressed from output
  # Also, reduce the output xtab starting column by one for every suppressed
  # colcol that comes before it
  outxtabcol = xtabcol
  for(i = 1; i <= length(colcolsarr); i++) {
    colcolshash[colcolsarr[i]] = 1
	if(colcolsarr[i] < xtabcol) --outxtabcol
  }
  # Similary for rowcols, but no reduction
  for(i = 1; i <= length(rowcolsarr); i++) {
    rowcolshash[rowcolsarr[i]] = 1
  }
  
  # First pass: read everything to get the column tuples
  outcol = outxtabcol
  j = 1
  if(header) getline line < FILENAME
  while(getline line < FILENAME) {
    split(line, row, FS)
	key = getKey(row, colcolsarr)
	if(!colkey2outcol[key]) {
	  colkeys[j++] = key
      colkey2outcol[key] = outcol++
	}
  }

  # Record the input col to output col mapping for extra ("more") columns
  adj = 0
  for(i = 1; i <= NF; i++) {
    # If it's a row-defining column, then it's not a more column, skip
    if(rowcolshash[i]) continue
	# If it's a col-defining column, then it's not a more column, and it
    # subtracts one from the adjustment
    else if(colcolshash[i]) {
      --adj
	  continue
	}
	# If it's the xtab col, then it's not a more column, and it affects the
	# adjustment
    if(i == xtabcol || colcolshash[i]) {
	  adj = adj - 1 + length(colkeys)
	  continue
	}
    morecolsin2out[i] = i + adj
  }

  # Print the header
  split("", outrow, "")
  j = 1;
  for(i = 1; i < xtabcol; i++) if(!colcolshash[i]) outrow[j++] = header ? $i : i
  for(i = 1; i <= length(colkeys); i++) outrow[j++] = colkeys[i]
  for(i = xtabcol + 1; i <= numcols; i++) if(!colcolshash[i]) outrow[j++] = header ? $i : i
  printRow()
  
  # Init prevrowkey
  prevrowkey = ""
  
  # If there was a header, then skip it before going on to the data
  if(header) getline
}
function initRow(rowkey, row) {
  prevrowkey = key
  
  split("", outrow, "")
  j = 1;
  for(i = 1; i < xtabcol; i++) if(!colcolshash[i]) outrow[j++] = row[i]
  for(i = 1; i <= length(colkeys); i++) outrow[j++] = ""
  for(i = xtabcol + 1; i <= numcols; i++) if(!colcolshash[i]) outrow[j++] = row[i]
}
function printRow() {
  printf("%s", outrow[1])
  for(i = 2; i <= length(outrow); i++) printf("\t%s", outrow[i])
  printf("\n")
}
function aggregate(old, new, fun) {
  if(old != "") {
    if(fun == "min") {
      if(old < new) new = old
	}
	else if(fun == "max") {
	  if(old > new) new = old
	}
	else if(fun == "sum") {
	  new += old
	}
	else if(fun == "first") {
	  new = old
    }
	else if(fun ~ /^cat$/) {
	  new = old new
	}
	else if(fun ~ /^cat./) {
	  new = old substr(fun, 4) new
	}
	# last/default: nothing to do
  }
  return new
}
# Now, for every row, including the first (this normal pass is the second pass)
{
  # Split into an array for getting keys
  split($0, row, FS)
  
  # Get row key, and see if it has changed
  rowkey = getKey(row, rowcolsarr)
  if(rowkey != prevrowkey) {
    if(prevrowkey) printRow()
	initRow(rowkey, row)
  }
  
  # Get colkey
  colkey = getKey(row, colcolsarr)
  
  # Determine new xtab value (default to value in current row, but update as
  # needed for aggregation)
  outcol = colkey2outcol[colkey]
  outrow[outcol] = aggregate(outrow[outcol], row[xtabcol], xtabfun)
  
  # Aggregate more columns
  for(incol in morecolsin2out) {
    outcol = morecolsin2out[incol]
	fun = morefunsarr[incol]
	if(!fun) fun = morefun
    outrow[outcol] = aggregate(outrow[outcol], row[incol], morefun)
  }
}
# In the end, print the final row
END {
  if(prevrowkey) printRow()
}
