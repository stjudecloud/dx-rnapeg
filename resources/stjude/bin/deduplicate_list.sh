#!/bin/bash
# De-duplicates the given list by removing duplicate elements and prints
# out the list containing unique elements to stdout.
#
# $1 = List from which duplicate elements are to be removed
# $2 = Delimiter used to separate elements (Ex: ",", ";", "|", etc.)

# Show usage information if insufficient parameters were sent
if [ "$#" -lt 2 ]; then about.sh $0; exit 1; fi

# Get parameter
list=$1
delimiter=$2

# Remove duplicates
echo $list | tr -s "$delimiter" "\n" | sort | uniq | tr "\n" "$delimiter" | sed 's/'$delimiter'$//'