#!/bin/gawk -f
# Reorders columns by sorting their values in the first row
#
# The output FS is set to the input FS, so if you leave it at default it will
# split on whitespace and output single space delimiters.  This is often not
# what you want to do, so typically, you will want to pass -F
#
# If you specify a negative right, then that is the number of columns NOT to
# sort at the right side (e.g. -v right=-1 leaves out the rightmost 1 column from
# reordering)
#
# Parameters:
# left: leftmost column to reorder (default 1)
# right: rightmost column to reorder (default: last column)
# sortargs: arguments to send to linux sort (optional, do not use -k)
# ordfile: file with column values in desired order (one per line, optional,
#   skips sort)
# defaultord: default ordinal to use if column header is not found in ordfile
#   (optional; default 99999 if ordfile, to sort to end)
BEGIN {
  OFS = FS
}
# The first row is where everything gets set up, since we read th values in the
# first row to determine the shuffling
NR == 1 {
  # Default in left and right if not specified
  if(!left) left = 1
  if(right <= 0) right += NF

  # Sort type to use (default of "" lets it be overridden by user in sortargs)
  sorttype = ""
  
  # Get mapping of new colnum to old colnum in new2old
  # If there is an ordfile, then use it
  if(ordfile) {
    ord = 1
    while(getline val < ordfile) val2ord[val] = ord++
	close(ordfile)
	if(defaultord = "") defaultord = 99999;
	# Like numeric, but sorts words to the end
	sorttype = "V"
  }
  # Using first line, determine the sort order
  cmd = "sort -k 1,1" sorttype " -k 2,2n"
  if(sortargs) cmd = cmd " " sortargs
  for(i = left; i <= right; i++) {
    ord = val2ord[$i]
	if(!ord) ord = defaultord $i
    printf("%s\t%d\n", ord, i) |& cmd
  }
  close(cmd, "to")
  new = left
  while((cmd |& getline line) > 0) {
    split(line, parts, "\t")
    new2old[new++] = parts[2]
  }
  close(cmd, "from")
  
  # Error checking
  if(new != right + 1) {
    print "Something we wrong, new = " new ", right = " right
	exit 1
  }
}
# Now, for every row, apply the column reordering
{
  for(i = left; i <= right; i++) saved[i] = $i
  for(i = left; i <= right; i++) $i = saved[new2old[i]]
  print
}
