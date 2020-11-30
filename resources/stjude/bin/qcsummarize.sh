#!/bin/sh
# Summarizes results for multiple test cases based on files containing the
# output.
#
# $1... = all arguments are file names

cat $* | awk '
  BEGIN { FS="\t"; OFS="\t"; cur=""; passed=0; failed=0; warned=0; }
  $1 == "HEADER" { cur=$2 ; desc=$3 }
  $1 == "SUMMARY" { summary=$3 }
  ($2 == cur) && ($1 ~ /(PASS)|(FAIL)|(WARN)/) {
    print $1, cur, desc, summary
  }
  ($2 == cur) && ($1 ~ /PASS/) { ++passed }
  ($2 == cur) && ($1 ~ /FAIL/) { ++failed }
  ($2 == cur) && ($1 ~ /WARN/) { ++warned }
  END {
    if(failed > 0) status = "FAIL";
    else if(warned > 0) status = "WARN";
    else status = "PASS";
    print status, "Summary", "Test cases: " (passed + failed + warned) "  Passed: " passed "  Failed: " failed "  Warnings: " warned;
  }'
