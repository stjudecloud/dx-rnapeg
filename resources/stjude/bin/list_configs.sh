#!/bin/bash
# Writes out all possible configurations for a category
#
# $1 = category name, or . for root category

dir=`which_config.sh $1`
for file in $dir/*.config.txt; do basename $file .config.txt; done
