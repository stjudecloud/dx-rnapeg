package LiftOverStrand;
# determine whether liftOver of a region is on the same strand

use strict;
use Carp qw(confess);

use Configurable;
use Exporter;
use FAI;
use TdtConfig;
use GenomeUtils qw(reverse_complement);
use MiscUtils qw(dump_die);
use List::Util qw(min);
use BLATer;

#use constant MIN_COMPARE_DISTANCE => 25;
#use constant MIN_COMPARE_DISTANCE => 60;
# need a good-length chunk particularly if mapping is via lower liftOver minMatch setting
#use constant MIN_TUPLES_TO_CALL => 20;

use constant MIN_COMPARE_DISTANCE => 25;
# tradeoff: want to evaluate region immediately around the site,
# OTOH need some space for comparison

use constant MIN_HIT_FRACTION => 0.75;
# what fraction of tuples is required to make a call

use constant MAX_MINOR_FRACTION => 0.20;
# what fraction of matches on the other strand is too many (e.g. reversible region)

my $TUPLE_SIZE = 12;
my $VERBOSE = 0;

my $BLAT_MIN_SCORE = 20;
my $BLAT_STEP_SIZE = 2;
# short reads
my $BLAT_MIN_IDENTICAL = 0.90;

@LiftOverStrand::ISA = qw(Configurable Exporter);
@LiftOverStrand::EXPORT_OK = qw();

use MethodMaker qw(
genome_from
genome_to
fais
chunk_cache_chromosomes
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->fais({});
  $self->chunk_cache_chromosomes(1);
  $self->configure(%options);
  return $self;
}

sub get_fai {
  my ($self, %options) = @_;
  my $genome = $options{"-genome"} || die;
  
  my $fais = $self->fais;
  my $fai = $fais->{$genome};
  unless ($fai) {
    my $config = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
    $fai = new FAI("-fasta" => $config->{FASTA} || die);
    $fai->chunk_cache_chromosomes(1) if $self->chunk_cache_chromosomes();
    $fais->{$genome} = $fai;
  }
  return $fai;
}

sub strand_check {
  my ($self, %options) = @_;

  my $is_same_strand;

  my $fai_from = $self->get_fai("-genome" => $self->genome_from);
  my $fai_to = $self->get_fai("-genome" => $self->genome_to);

  my $from_chr = $options{"-from-chr"} || die;
  my $from_start = $options{"-from-start"} || die;
  my $from_end = $options{"-from-end"} || die;

  my $to_chr = $options{"-to-chr"} || die;
  my $to_start = $options{"-to-start"};
  my $to_end = $options{"-to-end"};

  foreach ($to_start, $to_end) {
    die "to start/end undefined" unless defined $_;
  }

  my $dist1 = ($from_end - $from_start) + 1;
  my $dist2 = ($to_end - $to_start) + 1;

  if ($dist1 < MIN_COMPARE_DISTANCE or $dist2 < MIN_COMPARE_DISTANCE) {
    my $least = min($dist1, $dist2);
    my $needed = (MIN_COMPARE_DISTANCE - $least);
    my $flank = int($needed / 2) + 1;
    $from_start -= $flank;
    $from_end += $flank;
    $to_start -= $flank;
    $to_end += $flank;

    foreach ($from_start, $to_start) {
      $_ = 1 if $_ < 1;
    }
  }

  if ($fai_from->find_name($from_chr) and
      $fai_to->find_name($to_chr)) {

    my $chunk_from = $fai_from->get_chunk(
      "-id" => $from_chr,
      "-start" => $from_start,
      "-end" => $from_end,
	);

    my $chunk_to = $fai_to->get_chunk(
      "-id" => $to_chr,
      "-start" => $to_start,
      "-end" => $to_end,
	);
    foreach ($chunk_from, $chunk_to) {
      $_ = uc($_);
    }

    if ($chunk_from eq $chunk_to) {
      # simplest case: same strand
      $is_same_strand = 1;
    } else {
      my $chunk_from_rc = reverse_complement($chunk_from);
      if ($chunk_from_rc eq $chunk_to) {
	# simple case: opposite strand
	$is_same_strand = 0;
      } else {
	#
	# first lookup method: compare tuple counts.
	# if the target regions are of different lengths this method is required.
	#
	my $tuples_from = get_tuples(\$chunk_from);
	my $tuples_from_rc = get_tuples(\$chunk_from_rc);
	my $tuples_to = get_tuples(\$chunk_to);

	my $tuple_count = scalar keys %{$tuples_from};

	my $min_tuples_to_call = int($tuple_count * MIN_HIT_FRACTION);
	my $max_minor = int($tuple_count * MAX_MINOR_FRACTION);

	my $count_same = compare_tuples($tuples_from, $tuples_to);
	my $count_rc = compare_tuples($tuples_from_rc, $tuples_to);


	if ($VERBOSE) {
	  printf STDERR "from:%s from_rc:%s to:%s\n", $chunk_from, $chunk_from_rc, $chunk_to;
	  printf STDERR "to:%d-%d tuples:%d same:%d rc:%d\n", $to_start, $to_end, scalar(keys %{$tuples_from}), $count_same, $count_rc;
	}
	if ($count_same >= $min_tuples_to_call and
	    $count_same > $count_rc and
	    $count_rc <= $max_minor) {
	  $is_same_strand = 1;
	} elsif ($count_rc >= $min_tuples_to_call and
		 $count_rc > $count_same and
		 $count_same <= $max_minor) {
	  $is_same_strand = 0;
#	} else {
#	  $is_same_strand = "";
#	  dump_die(\%options, "insufficient data to decide", 1);
	} elsif (length($chunk_from) == length($chunk_to)) {
	  #
	  # tuples method fails: mismatches in target region?
	  # if the regions are of identical lengths, compare mismatch levels
	  #
	  my $mm_same = compare(\$chunk_from, \$chunk_to);
	  my $mm_rc = compare(\$chunk_from_rc, \$chunk_to);
	  my $len = length($chunk_from);
	  my $max_minor_mm = int($len * MAX_MINOR_FRACTION);

	  printf STDERR "string lookup: same=%d rc=%d\n", $mm_same, $mm_rc if $VERBOSE;
	  if ($mm_same != $mm_rc and
	      ($mm_same <= $max_minor_mm or $mm_rc <= $max_minor_mm)
	      # require relatively unambiguous result
	      ) {
	    if ($mm_same < $mm_rc) {
	      # fewer mismatches in + comparison
	      $is_same_strand = 1;
	    } else {
#	      dump_die(\%options, "hey now - $mm_same $mm_rc");
	    }
	  }
	}

	unless (defined $is_same_strand) {
	  #
	  #  all else fails: try blat
	  #  hg38->hg19 example: chr1.148972321.C.T
	  #
	  my $blat = new BLATer();
	  $blat->minScore($BLAT_MIN_SCORE);
	  $blat->stepSize($BLAT_STEP_SIZE);

	  my $parser = $blat->blat(
	    "-query" => {"genome_from" => $chunk_from},
	    "-database" => {"genome_to" => $chunk_to},
	      );
	  my $result = $parser->next_result();
	  # one result object per query sequence

	  if ($result) {
	    my @hsp;
	    while (my $hit = $result->next_hit()) {
	      # hits from this query to a database sequence
	      # (Bio::Search::Hit::HitI)
	      while (my $hsp = $hit->next_hsp) {
		push @hsp, $hsp if $hsp->frac_identical("query") >= $BLAT_MIN_IDENTICAL;
	      }
	    }
	    if (@hsp == 1) {
	      my $strand = $hsp[0]->strand();
	      if ($strand == -1) {
		$is_same_strand = 0;
	      } elsif ($strand == 1) {
		$is_same_strand = 1;
	      }
	    }
	  }
	  $is_same_strand = "" unless defined($is_same_strand);
	  printf STDERR "blat strand call for %s:%d-%d to %s:%d-%d: %s\n",
	  $from_chr, $from_start, $from_end, $to_chr, $to_start, $to_end, $is_same_strand;
	}

	die unless defined $is_same_strand;
#	dump_die(\%options, "no tuple hits to either") unless $count_same or $count_rc;
	printf STDERR "strand decision: %s\n", $is_same_strand if $VERBOSE;
      }
    }
  } else {
    $is_same_strand = "";
    # can't look up: mapped to a sequence not in target FASTA file
  }

  die unless defined $is_same_strand;
  return $is_same_strand;
}

sub get_tuples {
  # STATIC
  my ($string_ref) = @_;
  my $end = length $$string_ref;
  my $ei;
  my %tuples;
  for (my $si = 0; $si < $end; $si++) {
    # simpler way with regexp??
    $ei = $si + $TUPLE_SIZE - 1;
    last if $ei >= $end;
    my $tuple = substr($$string_ref, $si, $TUPLE_SIZE);
    $tuples{$tuple}++;
  }
#  foreach (keys %tuples) {
#    printf "%s: %d\n", $_, $tuples{$_};
#  }

  return \%tuples;
}

sub compare_tuples {
  # STATIC
  my ($tuples_from, $tuples_to) = @_;
  my $hits = 0;
  foreach my $t (keys %{$tuples_from}) {
    $hits++ if $tuples_to->{$t};
  }
  return $hits;
}


sub compare {
  my ($from, $to) = @_;
  die unless length $$from == length $$to;
  my $len = length $$from;
  my $mismatches = 0;
  for (my $i=0; $i < $len; $i++) {
    $mismatches++ unless substr($$from, $i, 1) eq substr($$to, $i, 1);
  }
  return $mismatches;
}

sub reference_swap_check {
  # 
  # check a single base in the two genomes for evidence of a swapped
  # reference/variant allele in the new genome.
  # FWIW: this likely does not interact well with strand-swapped situations!
  #
  my ($self, %options) = @_;
  my $from_chr = $options{"-from-chr"} || die;
  my $from_base = $options{"-from-pos"} || die;
  my $to_chr = $options{"-to-chr"} || die;
  my $to_base = $options{"-to-pos"} || die;
  my $alt_base = $options{"-alt-base"} || die;
  my $fai_from = $self->get_fai("-genome" => $self->genome_from);
  my $fai_to = $self->get_fai("-genome" => $self->genome_to);

  my $swapped;
  if ($fai_from->find_name($from_chr) and
      $fai_to->find_name($to_chr)) {
    my $chunk_from = $fai_from->get_chunk(
      "-id" => $from_chr,
      "-start" => $from_base,
      "-length" => 1
	);
    my $chunk_to = $fai_to->get_chunk(
      "-id" => $to_chr,
      "-start" => $to_base,
      "-length" => 1
	);

    foreach ($chunk_from, $chunk_to) {
      $_ = uc($_);
    }

    if ($chunk_from ne $chunk_to and
	$chunk_to eq $alt_base
	# base has changed between genomes, and new reference base is the old
	# variant base
	) {
      $swapped = 1;
    }
  }
  return $swapped;
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
