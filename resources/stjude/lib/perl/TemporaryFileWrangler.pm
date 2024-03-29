package TemporaryFileWrangler;
# create and clean up temporary files
# centralize local requirements, i.e. use of /bin/tmpnam
# MNE 8/2014

use strict;
use Carp qw(confess cluck);

use Configurable;
use Exporter;
use POSIX qw(tmpnam);

@TemporaryFileWrangler::ISA = qw(Configurable Exporter);
@TemporaryFileWrangler::EXPORT_OK = qw();

my %WARNED;

use MethodMaker qw(
use_mktemp
	auto_unlink
tempfiles
globfiles
verbose
mktemp_local_disk
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->tempfiles([]);
  $self->globfiles([]);
  $self->auto_unlink(1);
  $self->use_mktemp(1) unless $ENV{NO_MKTEMP};
  # SJ requirement, see
  # http://hc-wiki.stjude.org/display/compbio/Source+Code+Policies

  foreach my $dir ("/scratch_local",
		   "/scratch_space") {
    if (-d $dir and -w $dir) {
      $self->mktemp_local_disk($dir);
      last;
    }
  }
  # "mktemp --tmpdir=/scratch_local" will create a temporary directory under
  # /scratch_local. But it will need to be cleaned up by user..
  # - Zhou Mi, 11/1/2016
  #
  # "We just don't want users to fill up local disk on the head node by
  # accident. Can you use /scratch_space?"
  # - Zhou Mi, 1/27/2017
  $self->configure(%options);
  return $self;
}

sub get_tempfile {
  # TO DO:
  my ($self, %options) = @_;
  my $append = $options{"-append"};
  # optionally append a string
  my $fn;

  if ($self->use_mktemp) {
    my $cmd = "mktemp";
    my $local_dir = $self->mktemp_local_disk();
    if ($local_dir) {
      if (-d $local_dir) {
	if (-w $local_dir) {
	  $cmd .= sprintf ' --tmpdir=%s', $local_dir;
	} else {
	  printf STDERR "WARNING: specified temp root %s is not writable\n", $local_dir;
	}
      } else {
	printf STDERR "WARNING: specified temp root %s does not exist!\n", $local_dir unless $WARNED{$local_dir};
	$WARNED{$local_dir} = 1;
      }
    }
    $fn = `$cmd`;
    chomp $fn;
    die sprintf "can't create tempfile, file=%s exit=%s", $fn || "null", $? if $? or not($fn);
#    cluck "create tempfile $fn";
  } else {
    $fn = tmpnam();
  }

  $self->add_tempfile($fn);
  # a zero-length version of this file always seems to be created
  # even if user only wants the appended version below
  $self->add_glob($fn) if $options{"-glob"};
  # delete all files starting with this name
  # (e.g. temporary reference sequence file that then has index files added)

  if ($append) {
    # add user-provided tag
    $fn .= $append;
    $self->add_tempfile($fn);
  }

  if (-s $fn) {
    printf STDERR "WARNING: tempfile %s already exists and is size > 0!\n", $fn;
    unlink $fn;
    # TO DO: maybe just skip and call tempfile code again? (w/max retries?)
    confess "WTF: temp outfile $fn already exists and is size > 0" if -f $fn;
  }

#  cluck($fn);

  return $fn;
}

sub add_tempfile {
  my ($self, $fn) = @_;
  push @{$self->tempfiles}, $fn;
}

sub add_glob {
  my ($self, $fn) = @_;
  push @{$self->globfiles}, $fn;
}

sub DESTROY {
  # auto-clean tempfiles
  my ($self) = @_;
  my $tempfiles = $self->tempfiles();
  my $globfiles = $self->globfiles();
  my @tmp;
  push @tmp, @{$tempfiles};
  foreach my $g (@{$globfiles}) {
    push @tmp, glob($g . "*");
  }

  if (@tmp) {
    if ($self->auto_unlink()) {
      printf STDERR "TemporaryFileWrangler DESTROY: deleting %s\n", join " ", @tmp if $self->verbose;
      unlink @tmp;
    } else {
      printf STDERR "\nTemporaryFileWrangler.pm: NOT unlinking tempfiles %s\n", join ",", @tmp;
    }
  }
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
