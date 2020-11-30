#!/bin/env perl
# submit jobs for a list of files, with awareness of cluster status/completion

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use FileUtils qw(read_simple_file);
use MiscUtils qw(dump_die);
use Cluster;
use File::Basename;

my %FLAGS;
GetOptions(\%FLAGS,
	   "-file-glob=s",
	   # specify input files by glob
	   "-files=s",
	   "-cmds=s",

	   "-ram=i",
	   "-template=s",

	   "-out-suffix=s",
	   "-debug",
	  );

my $ram = $FLAGS{ram} || die "-ram";
my $template = $FLAGS{template};
if ($template) {
  die "template needs %s" unless $template =~ /%s/;
}
my $cmds = $FLAGS{cmds};
die "need -template or -cmds" unless $template or $cmds;

my $out_suffix = $FLAGS{"out-suffix"};
unless ($cmds) {
  die "-out-suffix" unless $out_suffix;
  die "no leading period in -out-suffix" if $out_suffix =~ /^\./;
}

my @infiles;
if (my $g = $FLAGS{"file-glob"}) {
  @infiles = glob($g);
} elsif (my $lf = $FLAGS{files}) {
  my $set = read_simple_file($lf);
  @infiles = @{$set};
} elsif ($cmds) {
  my $set = read_simple_file($cmds);
  @infiles = @{$set};
} else {
  die "specify -files LISTFILE -file-glob [pattern]\n";
}

die "no files" unless @infiles;

my $counter=0;
foreach my $infile (@infiles) {
  $counter++;
  my $cmd;
  my $outfile;
  if ($cmds) {
    $cmd = $infile;
    $outfile = sprintf "bogus_%d.out", ++$counter;
  } else {
    $cmd = sprintf $template, $infile;
    $outfile = basename($infile) . "." . $out_suffix;
  }

  my $c = new Cluster(
		      "-memory" => $ram,
		      "-command" => $cmd,
		      "-project" => "PCGP",
		      "-outfile" => $outfile,

		      "-debug" => $FLAGS{debug},
		     );
  $c->run();

}

