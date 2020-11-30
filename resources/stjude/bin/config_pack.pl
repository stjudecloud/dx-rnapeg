#!/bin/env perl
# pack up genome config data for an application, e.g. for cloud use
# MNE 12/2016

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version
use Carp qw(confess);
use Getopt::Long;
use File::Basename;
use File::Path;

use MiscUtils qw(log_message);
use FileUtils qw(read_simple_file);
use DelimitedFile;
use Reporter;
use TdtConfig;

use constant WC_ROOT => "/research/rgs01/resgen/dev/wc/";
use constant CLOUD_ROOT => "/usr/local/tartan/index/reference/Homo_sapiens/";

#my $COMPRESSION_LEVEL = 3;
my $COMPRESSION_LEVEL = 2;
# - gzip -9 is VERY slow
# - level 2 or 3 is still pretty fast but saves space compared to 1
# - many large resources are compressed already (e.g. tabix files),
#   so little point in 3+

my %FLAGS;
my @clopts = (
	      "-genome=s",
	      # e.g. GRCh37-lite
	      "-package=s",
	      # compbio research cluster code package, e.g. snv-annovar
	      "-wc=s",
	      # manual override only

	      "-code=s",
	      # rather than packing up config data, put module code
	      # in cloud resources dir of specified app

	      "-level=i" => \$COMPRESSION_LEVEL,
	      "-bzip2",
	      "-no-archive",
	     );
GetOptions(
	   \%FLAGS,
	   @clopts
	  );

my $package = $FLAGS{package} || die "-package";

my $tr = $ENV{TARTAN_ROOT} || die "TARTAN_ROOT not set, run setcbenv";
printf STDERR "TARTAn root: %s\n", $tr;
$tr =~ /\/(\w+)\/tartan/ || die;
my $tartan_env = $1;

my $wc_dir = $FLAGS{wc};
unless ($wc_dir) {
#  my $uid = getlogin || getpwuid($<) || die "can't determine user id";
  my $uid = getpwuid($<) || die "can't determine user id";
  # getlogin() returns "root" on compute nodes (?)
  $wc_dir = sprintf '%s/%s', WC_ROOT, $uid;
}
die "where is $wc_dir" unless -d $wc_dir;

if ($FLAGS{code}) {
  pack_code();
} else {
  pack_config_data();
}

sub pack_config_data {
  my $genome = $FLAGS{genome} || die "-genome";
  my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";

  my $wc_config = sprintf '%s/cluster_code/trunk/configs/data/genome/%s.config.txt', $wc_dir, $genome;
  die "where is $wc_config" unless -s $wc_config;
  # working copy of genome config contains human-organized paths

  my $f_config_variables = sprintf '%s/cluster_code/trunk/%s/genome_files.txt', $wc_dir, $package;
  die "where is $f_config_variables" unless -s $f_config_variables;

  my $v_wanted = read_simple_file($f_config_variables, "-hash1" => 1);
  # config variables touched by the application to pack up

  open(IN, $wc_config) || die;
  my $cloud_root = sprintf '%s/%s', CLOUD_ROOT, $genome;
  $cloud_root =~ s/\/\//\//g;

  my $pack_subdir = "export";
  rmtree($pack_subdir);
  die if -d $pack_subdir;

  my $f_configs = sprintf '%s/usr/local/configs/data/genome/%s', $pack_subdir, basename($wc_config);
  # data/genome: hack, need cleaner way
  mkpath(dirname($f_configs));
  my $wf = new WorkingFile($f_configs);
  my $fh = $wf->output_filehandle();

  #
  #  create links to touched resources and build config file:
  #
  while (my $line = <IN>) {
    chomp $line;
    next if $line =~ /^#/;
    next unless $line =~ /\w/;
    my @f = split /\t/, $line;
    printf STDERR "WARNING: malformed line with %d fields: %s\n", scalar(@f), $line unless @f == 2;
    my ($variable, $value) = @f;
    if ($v_wanted->{$variable}) {
      printf STDERR "processing %s\n", $variable;
      my $rel_path = $value;
      if ($rel_path =~ /^\$/) {
	$rel_path =~ s/^\$\w+\/// || die "can't strip variable from $value";
	# get human-organized relative path
	$rel_path =~ s/\/$//;
	# remove trailing slash for directories, interferes w/symlinking
      }

      my $deployed_value = $config_genome->{$variable} || die;
      my $cloud_value;
      my $need_link;
      my $is_dir;

      if (-f $deployed_value or -d $deployed_value) {
	$cloud_value = sprintf '%s/%s', $cloud_root, $rel_path;
	$need_link = 1;
	#      $is_dir = 1 if -d $deployed_value;
      } else {
	# string
	$cloud_value = $deployed_value;
      }

      printf $fh "%s\n", join "\t", $variable, $cloud_value;

      if ($need_link) {
	my $local = $pack_subdir . $cloud_value;

	unlink $local;
	mkpath(dirname($local));
	#      printf STDERR "symlink %s => %s\n", $deployed_value, $local;
	symlink($deployed_value, $local) || die "can't symlink $deployed_value => $local: $!";

	if ($local =~ /\.fa$/i) {
	  # special handling for FASTA: also include samtools .fai index
	  my $fai = $deployed_value . ".fai";
	  if (-s $fai) {
	    unlink $local . ".fai";
	    symlink($deployed_value . ".fai", $local . ".fai") || die;
	  }
	}
      }
    }
  }
  $wf->finish();

  if ($FLAGS{"no-archive"}) {
    printf STDERR "skipping archive\n";
  } else {
    #
    #  build archive:
    #
    my ($compressor, $suffix);
    if ($FLAGS{bzip2}) {
      $compressor = "bzip2";
      $suffix = "bz2";
    } else {
      $compressor = "gzip";
      $suffix = "gz";
    }

    my $f_archive = sprintf "bundle_%s_%s_%s.tar.%s", $tartan_env, $package, $genome, $suffix;

    my $cmd = sprintf '(cd %s/; tar -hvcf - . |%s -%d) > %s', $pack_subdir, $compressor, $COMPRESSION_LEVEL, $f_archive;
    log_message("building archive: $cmd");
    system $cmd;
    die "$cmd exited with $?" if $?;
    log_message("done");
  }
}

sub pack_code {
  my $cloud_pkg = $FLAGS{code} || die;

  my $cloud_pkg_dir = sprintf '%s/apps/dnanexus/%s', $wc_dir, $cloud_pkg;
  die "where is $cloud_pkg_dir" unless -d $cloud_pkg_dir;

  my $cloud_resources_target = sprintf '%s/resources/stjude/', $cloud_pkg_dir;
  mkpath($cloud_resources_target);

  my $code_wc = sprintf '%s/cluster_code/trunk/%s', $wc_dir, $package;
  die "where is $code_wc" unless -d $code_wc;
  chdir($code_wc);

  my $cmd = sprintf 'sjcb_build.sh --and-deps --package release:%s', $cloud_resources_target;

  system $cmd;
  die "$cmd failed with $?" if $?;

  # TO DO:
  # remove emacs backup files
}
