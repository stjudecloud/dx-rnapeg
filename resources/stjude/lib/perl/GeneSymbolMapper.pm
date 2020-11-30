package GeneSymbolMapper;
# attempt to translate gene symbols in one set to another
# MNE 4/2014

# tests:
#
# SV:
#
#
# CNV analysis using SV genes in GENE_EXON_REGION:
# MADH4: should be SMAD4 (later symbol)
# WTX:  updated symbol is AMER1 but GENE_EXON_REGION uses FAM123B (older)
#

use strict;
use Carp qw(confess);

use Configurable;
use Exporter;

use MiscUtils qw(dump_die);
use HGNCParser;
use GeneSymbolStandardizer;

@GeneSymbolMapper::ISA = qw(Configurable Exporter);
@GeneSymbolMapper::EXPORT_OK = qw();

use MethodMaker qw(
	genes
	hgnc
	gss

	hgnc_file
	eg_file
	refseq2gene

problem
approved_symbol

hgnc_synonym_prune
hgnc_synonym_enable
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->genes({});
  $self->hgnc_synonym_enable(1);
  $self->hgnc_synonym_prune(1);
  # enable HGNC synonyms by default, but prune some synonyms
  # to reduce risk of dubious matches.
  #
  # example:
  # - existing germline reviewable list contains FANCA
  # - new reviewable gene FAH
  # - both FAH and FANCA are approved symbols
  # - however, FAH is also be a synonym for FANCA
  # - these are two completely different genes, a synonym match
  #   is not strong enough evidence to reject gene as a duplicate
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self) = @_;
  $self->refseq2gene({});

  my $hgnc = new HGNCParser("-file" => $self->hgnc_file || die "-hgnc_file");
  $hgnc->prune_synonyms() if $self->hgnc_synonym_prune();
  $self->hgnc($hgnc);
  # HGNC/HUGO

  $self->gss(new GeneSymbolStandardizer("-filename" => $self->eg_file || die "-eg_file"));
  # Entrez Gene
}

sub add_gene {
  my ($self, %options) = @_;
  my $gene = $options{"-gene"} || confess "-gene";
  my $refseq = $options{"-refseq"};
  # optional
  my %record;
  $record{gene} = $gene;
  $record{refseq} = $refseq;

  if ($refseq) {
    my $r2g = $self->refseq2gene();
    if (my $existing = $r2g->{$refseq}) {
      printf STDERR "WARNING: duplicate refGene entry for %s: old:%s new:%s\n", $refseq, $existing, $gene unless $existing eq $gene;
      # can happen for hack NG_ entries, see COMPBIO-2772
    } else {
      $r2g->{$refseq} = $gene;
    }
  }

  push @{$self->genes->{$gene}}, \%record;
}


sub reset {
  my ($self) = @_;
  $self->problem("");
  $self->approved_symbol("");
  # reset before sift_hgnc_hits() called for the first time
}

sub resolve {
  my ($self, %options) = @_;
  my $sym = $options{"-symbol"} || die "-symbol";
  my $result;
  $self->reset();

  if (0 and $self->contains($sym)) {
    # target is the same in both lists.
    # normally we'd check before calling this routine and this would
    # be a fatal error as a lookup is not needed, and we don't
    # want to require database matching/filtering if so.
    #
    # However for some use cases (e.g. universal somatic gene list
    # construction) it's desirable to run these symbols through
    # to determine the HUGO symbol for each, i.e. to detect
    # synonymous symbols in combined list.
    die "same";
  }

  my $hgnc = $self->hgnc();
  my $problem;

  if (0) {
    my $test = $hgnc->find("-symbol" => $sym,
			   "-approved" => 1);
    dump_die($test->[0]);
  } 

  ($result, $problem) = $self->sift_hgnc_hits(%options,
					      "-hits" =>
					      $hgnc->find("-symbol" => $sym,
							  "-approved" => 1)
      ) unless $result or $problem;
  # ACKR3: a SV gene not found in GENE_EXON_REGION,
  # but NM_020311 and older symbol CXCR7 are

  ($result, $problem) = $self->sift_hgnc_hits(%options,
					      "-hits" =>
					      $hgnc->find("-symbol" => $sym,
							  "-previous" => 1)
      ) unless $result or $problem;
  # CEP1 to GENE_EXON_REGION via NM_007018:
  # CEP1 -> CNTRL -> NM_007018 -> CEP110
  # this one is tricky because CEP1 is also a synonym for CDC42EP1.
  # avoid by using previous symbol lookups before synonym lookups.

  if ($self->hgnc_synonym_enable) {
    ($result, $problem) = $self->sift_hgnc_hits(%options,
						"-hits" =>
						$hgnc->find("-symbol" => $sym,
							    "-synonym" => 1)
	) unless $result or $problem;
  }

  unless ($result or $problem) {
    # try Entrez Gene if all else fails, occasionally has
    # an entry (e.g. ALO17 -> RNF213)
    my $syns = $self->gss->find($sym);
    if ($syns and ref $syns) {
      if (@{$syns} == 1 and $self->contains($syns->[0])) {
	($result) = @{$syns};
      } else {
	$problem = 1;
      }
    }
  }

  printf STDERR "no luck for %s\n", $sym unless $result;


#     if ($hugo_hits) {
#       if (@{$hugo_hits} == 1) {
# 	# only 1 result
# 	my $approved = $hugo_hits->[0]->{"Approved Symbol"} || die;
# 	if ($self->contains($approved)){
# 	  # given symbol translates to an approved symbol which
# 	  # is found in the match list
# 	  $result = $approved;
# 	} else {
# 	  die "can't find approved $approved";
# 	  # try older symbols??
# 	}
#       } else {
# 	# ambiguous matches, e.g.
# 	# HSPCA, which hits HSP90AA1 and HSP90AA2.
# 	my @try = map {$_->{"Approved Symbol"}} @{$hugo_hits};
# #	die @try;
# 	my @matches;
# 	foreach my $g (@try) {
# 	  push @matches, $g if $self->contains($g);
# 	}

# 	if (@matches == 1) {
# 	  ($result) = @matches;
# 	} else {
# 	  die "not matches or ambiguous";
# 	}
#       }
#     }

  return $result;
}

sub contains {
  my ($self, $sym) = @_;
  return $self->genes->{$sym} ? 1 : undef;
}

sub sift_hgnc_hits {
  my ($self, %options) = @_;
  die unless exists $options{"-hits"};
  my $hugo_hits = $options{"-hits"};
  my $result;
  my %found;
  if ($hugo_hits) {

    my %approved;

    foreach my $hh (@{$hugo_hits}) {
      my $approved = $hh->{"Approved Symbol"} || die;
      $approved{$approved} = 1;

      my $match;
      if ($self->contains($approved)) {
	# latest approved symbol is in the target set
	$found{$approved} = 1;
	$match = 1;
      }

      if (not($match) and $hh->{refseqs_hash}) {
	# the latest approved symbol isn't in the target set,
	# however it does contain NM records, try that lookup.
	foreach my $rs (keys %{$hh->{refseqs_hash}}) {
	  my $hit = $self->refseq2gene->{$rs};
	  if ($hit) {
	    $found{$hit} = 1;
	    $match = 1;
	  }
	}
      }

      if (not($match) and my $ps = $hh->{"Previous Symbols"}) {
	# see if any of the previous symbols for this entry 
	# are found in the target set, e.g.
	#   1. source symbol (SV config): CMKOR1 ->
	#   2. approved symbol (FB): ACKR3 ->
	#   3. target symbol (GER): CXCR7
	my @prev = split /,\s*/, $ps;
	foreach my $sym (@prev) {
	  if ($self->contains($sym)) {
	    $found{$sym} = 1;
	    $match = 1;
	  }
	}
      }
    }

    if (%approved and !$self->approved_symbol()) {
      my $list = join ",", sort keys %approved;
      $self->approved_symbol($list);
      # only record first/best match
      # there may be more than one, e.g. lookup of SIL:
      # - match found via synonyms only
      # - hits multiple genes: STIL, PMEL
    }

  }

  my $problem;

  if (my @found = keys %found) {
    if (@found == 1) {
      ($result) = keys %found;
    } else {
      $problem = "problematic_ambiguity";
      printf STDERR "ERROR: %s for %s\n", $problem, $options{"-symbol"};
    }
  }

  return ($result, $problem);
}

sub populate_refflat {
  my ($self, %options) = @_;
  my $rf = $options{"-refflat"} || die "-refflat";

  open(RFTMP, $rf) || die;
  my $genes = 0;
  while (<RFTMP>) {
    chomp;
    my @f = split /\t/, $_;
    my $gene = clean_sharp_gene_symbol($f[0]);
    if ($gene) {
      $genes++;
      my $nm = $f[1];
      $self->add_gene(
		     "-gene" => $gene,
		     "-refseq" => $nm
		    );
    } else {
#      printf STDERR "no gene symbol in refFlat entry %s\n", $_;
      # e.g. _locPar NR_001526-locPar
    }
  }
  close RFTMP;
  die "no genes" unless $genes;
}


sub clean_sharp_gene_symbol {
  # in Michael Rusch's "sharp" versions of refFlat,
  # ambiguous loci assigned "_loc" suffix, e.g. DUX4_locF.
  # These will also be reported by FusionBuilder with this suffix,
  # so need to remove it to get the raw gene symbol.
  my ($gene) = @_;
  if ($gene =~ /_loc/) {
    my $before = $gene;
    $gene =~ s/_loc.*$//;
#    printf STDERR "stripping %s => %s\n", $before, $gene;
  }
  return $gene;
}

sub populate_gene_exon_region {
  #
  # symbols reported by GENE_EXON_REGION
  #
  my ($self, %options) = @_;
  my $ger_dir = $options{"-dir"} || die;
  my @files = glob("$ger_dir/*txt");
  die unless @files;
  foreach my $f (@files) {
    open(GERTMP, $f) || die;
    while (<GERTMP>) {
      my @f = split /\t/, $_;
      my ($gene, $nm) = @f[0,2];
      $self->add_gene(
		     "-gene" => $gene,
		     "-refseq" => $nm
		    );
    }
  }
  close GERTMP;
}

sub populate_chr2gene {
  my ($self, %options) = @_;
  my $dir = $options{"-dir"} || die;
  my @files = glob($dir . "/*gene.txt");
  die "no files in $dir" unless @files;
  my %saw;
  foreach my $f (@files) {
    open(CGTMP, $f) || die;
    while(<CGTMP>) {
      chomp;
      my @f = split /\t/, $_;
      die unless @f == 4;
      my @list = split /\|/, $f[0];
      my ($gene, $nm) = @list;
      unless ($saw{$gene}) {
	$self->add_gene(
			"-gene" => $gene,
			"-refseq" => $nm
		       );
	$saw{$gene} = 1;
      }
    }
  }
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
