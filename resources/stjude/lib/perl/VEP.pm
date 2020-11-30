package VEP;
# wrapper/prep for VEP (Variant Effect Predictor) 
# http://www.ensembl.org/info/docs/tools/vep/script/index.html

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use TemporaryFileWrangler;
use FileUtils qw(find_binary);
use WorkingFile;

@VEP::ISA = qw(Configurable Exporter);
@VEP::EXPORT_OK = qw();

use MethodMaker qw(
	variants
cache_dir
fasta
v2vep
tfw
vep_in
vep_out
vep_command
prep_only
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->variants([]);
  $self->tfw(new TemporaryFileWrangler());
  $self->configure(%options);
  return $self;
}

sub add_variant {
  my ($self, $v) = @_;
  push @{$self->variants}, $v;
}

sub run_vep {
  my ($self, %options) = @_;

  find_binary("variant_effect_predictor.pl", "-die" => 1);
  my $tfw = $self->tfw;
  my $vep_in = $tfw->get_tempfile("-append" => ".vep");
  my $vep_out = $vep_in . ".out";
  $tfw->add_tempfile($vep_out);

  $self->vep_in($vep_in);
  $self->vep_out($vep_out);

  my $wf = new WorkingFile($vep_in);
  my $fh = $wf->output_filehandle;

  my %key2v;
  #
  #  prep input file:
  #
  # http://www.ensembl.org/info/docs/tools/vep/vep_formats.html#input
  foreach my $v (@{$self->variants}) {
    my $key = $v->get_snv4();
    die "duplicate input $key" if $key2v{$key};
    # sanity check: shouldn't happen
    $key2v{$key} = $v;
    my $chr = $v->reference_name;
    my $ref_base = $v->reference_allele;
    my $var_base = $v->variant_allele;
    my $strand = "+";
    my $id = $key;
    my ($start, $end);
    if ($v->is_substitution() or $v->is_deletion() or $v->is_complex()) {
      $start = $v->start;
      $end = $v->end;
    } elsif ($v->is_insertion()) {
      # "An insertion (of any size) is indicated by start coordinate = end coordinate + 1."
      $start = $v->start + 1;
      $end = $v->end;
    } else {
      die "unhandled variant type";
    }

    printf $fh "%s\t%d\t%d\t%s/%s\t%s\t%s\n",
    $chr,
    $start,
    $end,
    $ref_base,
    $var_base,
    $strand,
    $id;
  }
  $wf->finish();

  return if $self->prep_only();

#  printf "%s\n", $vep_in; sleep 60;

  my $fasta = $self->fasta() || die "-fasta";
  my $cache_dir = $self->cache_dir() || die "-cache_dir";

  my $cmd = sprintf 'variant_effect_predictor.pl --quiet --no_stats --refseq --hgvs --force_overwrite --offline --fasta %s -i %s -o %s --dir %s --dir_cache %s --dir_plugins %s',
  $fasta,
  $vep_in,
  $vep_out,
  $cache_dir,
  $cache_dir,
  $cache_dir;
  
  $self->vep_command($cmd);

  system($cmd);
  die "$cmd exited with $?" if $?;

  #
  #  parse output:
  #
  my %v2vep;

  if (-s $vep_out) {
    open(VEPTMP, $vep_out) || die;
    my @headers;
    while (<VEPTMP>) {
      chomp;
      if (/^##/) {
	next;
      } elsif (/^#/) {
	die if @headers;
	s/^#//;
	@headers = split /\s+/, $_;
      } else {
	my @f = split /\s+/, $_;
	die "row/header mismatch" unless @f == @headers;
	my %r;
	@r{@headers} = @f;

	my $key = $r{Uploaded_variation} || die;
	my $v = $key2v{$key} || die;
	push @{$v2vep{$v}}, \%r;
      }
    }
  } else {
    die "no outfile $vep_out" unless -s $vep_out;
  }

  $self->v2vep(\%v2vep);
}

sub get_results {
  my ($self, $v) = @_;
  # $v = same Variant.pm reference as added
  return $self->v2vep->{$v};
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
