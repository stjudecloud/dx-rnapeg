#!/usr/bin/env perl
# multiplex processing of a tab-delimited file
# MNE 4/2015

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use File::Basename;
use POSIX qw(ceil);

use MiscUtils qw(dump_die build_argv_list get_core_count);
use DelimitedMux;
# use DelimitedFile;
# use Reporter;

my %FLAGS;
GetOptions(\%FLAGS,
	   "-file=s",
	   "-files=s",
	   "-glob=s",

	   "-template=s",
	   "-suffix=s",
	   "-ram=i",
	   "-wait=i",

	   "-count=i",
	   "-jobs=i",
	   "-percent=s",
	   "-pool=i",

	   "-clean=s",
	   "-glob=s",

	   "-parallel",
	   # use GNU parallel rather than submitting jobs to cluster

	   "-auto-core",
	   "-force-cores=i",
	   "-force-lines=i",
	   "-out=s",
	  );

my $files = build_argv_list("-flags" => \%FLAGS,
			    "-single" => "file",
			    "-set" => "files",
			    "-glob" => "glob",
			   );

my $RAM = $FLAGS{ram};
die "-ram" if $FLAGS{parallel} ? 0 : not($RAM);
my $template = $FLAGS{template} || die "-template";

#my $suffix = "annovar_merged.tab";
my $suffix = $FLAGS{suffix} || die "-suffix";
$suffix =~ s/^\.//;

my $outfile_manual = $FLAGS{out};
die "-out requires single input file" if $outfile_manual and @{$files} != 1;

foreach my $f (@{$files}) {
  my $of = $outfile_manual || basename($f) . "." . $suffix;
  my $count;
  if ($FLAGS{"auto-core"}) {
    my $preferred = $FLAGS{count} || die "-auto-core requires -count PREFERRED_ROW_COUNT";
    my $core_count = $FLAGS{"force-cores"} || get_core_count();
    my $line_count = $FLAGS{"force-lines"} || count_data_lines($f);

    my $job_count = ceil($line_count / $preferred);

    if ($job_count < $core_count) {
      # using the preferred count will result in fewer jobs than
      # available cores.  Use a smaller count to keep all cores busy.
      $count = ceil($line_count / $core_count);
    } else {
      $count = $preferred;
    }

    $job_count = ceil($line_count / $count);
    if ($job_count - $core_count <= 1) {
      # if specified count results in one more job than available cores,
      # increase slightly count slightly to spread load more evenly.
      # helpful esp. for jobs w/heavy startup overhead (e.g. medal ceremony)
      $count = ceil($line_count / $core_count);
    }

    $job_count = ceil($line_count / $count);
    printf STDERR "auto-core: %d jobs of %d lines\n", $job_count, $count;
  } elsif ($count = $FLAGS{count}) {
    # manually specified
  } elsif (my $jobs = $FLAGS{jobs} or my $percent = $FLAGS{percent}) {

    my $lines = count_data_lines($f);

    if ($jobs) {
      $count = int($lines / $jobs) + 1;
    } else {
      $count = 
      die "percent not yet implemented"
    }
  } else {
    die "specify -count|-jobs|-percent";
  }

  if (-s $of) {
    printf STDERR "%s already exists\n", $of;
  } else {
    my $dm = new DelimitedMux();
    if (my $t = $FLAGS{clean}) {
      $dm->cleanup_type($t);
    }

    $dm->split_file(
		    "-file" => $f,
		    "-lines" => $count,
		   );

    $dm->run_jobs(
		  "-template" => $template,
		  "-ram" => $RAM,
		  "-out-suffix" => $suffix,
		  "-pool" => $FLAGS{pool},
		  "-wait" => $FLAGS{wait},
		  "-parallel" => $FLAGS{parallel},
		 );

    $dm->join_files("-out" => $of);
  }
}

sub count_data_lines {
  my ($f) = @_;
  open(TMP, $f) || die "can't open $f";
  my $lines = 0;
  while (<TMP>) {
    $lines++;
  }
  $lines--;
  # -1 for required header line
  return $lines;
}
