#!/usr/bin/env perl
# annotate raw junction count files with gene name
# this is a stripped-down version of junction2gene.pl that:
# - avoids using the DBI.pm module which seems to have problems with
#   Docker on DNAnexus
# - only supports junction files, not .bed files

use strict;
use warnings;
use 5.10.1;
# helpful hint if attempted w/older version

use Getopt::Long;

use DelimitedFile;
use FileUtils qw(read_simple_file);
use GeneAnnotationFlat;

my %FLAGS;
my $FORCE = 1;
GetOptions(\%FLAGS,
	   "-junctions=s",
	   # single junction/count file
	   "-list=s",
	   # list of files

	   "-glob",
	   # glob for junction files
	   "-glob-pattern=s",

#	   "-keep-unknown",
	   "-genome=s",
	   "-no-config",

	   "-preserve",
	   # preserve existing annotations for junctions already
	   # annotated as known (just add gene symbols for novel)

	   "-gene=s",
	   # restrict processing to specified gene symbol

	   "-force=i" => \$FORCE,

	   "-reference-transcript=s",

#	   "-extract-excerpt",
	   # create an excerpt file from existing already-annotated
	   # output file

	   "-gene-exon-region-dir=s",
	   # annotate from GENE_EXON_REGION directory

	   "-refgene=s",
	   # annotate from refGene.txt

	   "-strand=s",

	  );

my $GA;
my $genome = $FLAGS{genome} || die "specify -genome";
# required for chromosome disambiguation code, e.g.
# in human 23 = X, but for zebrafish 23 = 23
my $strand_filter = $FLAGS{strand};

if (my $ga_dir = $FLAGS{"gene-exon-region-dir"}) {
  $GA = new GeneAnnotationFlat(
			   "-style" => "gene_exon_region",
			   "-gene_exon_region_dir" => $ga_dir,
			   "-ignore_non_coding" => 0,
			   "-genome" => $genome,
			   "-strand" => $strand_filter
			  );
} elsif (my $rg = $FLAGS{"refgene"}) {
  $GA = new GeneAnnotationFlat(
			   "-style" => "refgene_flatfile",
			   "-refgene_flatfile" => $rg,
			   "-ignore_non_coding" => 0,
			   "-genome" => $genome,
			   "-strand" => $strand_filter
			  );
} else {
  die "no DBI options in this script, see junction2gene.pl";
}

#
#  map junctions to genes:
#
my @junction_files;
if (my $fn = $FLAGS{junctions}) {
  push @junction_files, $fn;
} elsif (my $l = $FLAGS{list}) {
  my $list = read_simple_file($l);
  push @junction_files, @{$list};
} elsif ($FLAGS{glob}) {
  @junction_files = glob("*.junctions.tab");
  die "no matching files" unless @junction_files;
} elsif (my $p = $FLAGS{"glob-pattern"}) {
  @junction_files = glob($p);
  die "no matching files" unless @junction_files;
} else {
  die "specify -junctions or -glob\n";
}
die unless @junction_files;

foreach my $jf (@junction_files) {
  process_junction_file($jf);
}

sub process_junction_file {
  my ($junction_file) = @_;

  my $restrict_gene = $FLAGS{"gene"};
  my $restrict_ref = $FLAGS{"reference-transcript"};

#  my $suffix = $restrict_gene ? "annotated." . $restrict_gene . ".tab" : "annotated.tab";
  my @suffix = "annotated";
  push @suffix, $restrict_gene if $restrict_gene;
  push @suffix, $restrict_ref if $restrict_ref;
  push @suffix, "tab";

  my $suffix = join ".", @suffix;


#  my $outfile = $junction_file . ".annotated.tab";
  my $outfile = $junction_file . "." . $suffix;

  return if -s $outfile and !$FORCE;

  my $df = new DelimitedFile(
			     "-file" => $junction_file,
			     "-headers" => 0,
			    );
  $df->delimiter("\t");

  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   junction
					   count
					   genes
					)
				      ]
			);


  my $no_gene = 0;
  my $multi_genes = 0;
  my $known_mode = 0;
  my $preserve_mode = 0;
  my @headers;
  my %ignore_warn;
  my $transcript_index = -1;
  while (my @row = $df->next()) {
    my ($junction, $count) = @row;
    if ($junction eq "junction") {
      # header line (from SpliceReadReport -annotate mode)
      $known_mode = 1;
      if ($FLAGS{preserve}) {
	$preserve_mode = 1;
	@headers = @row;
	if ($restrict_ref) {
	  for (my $i = 0; $i < @headers; $i++) {
	    $transcript_index = $i if $headers[$i] eq "transcripts";
	  }
	}
      } else {
	@headers = qw(
		      junction
		      count
		      known
		      genes
		   );
      }
      $rpt->labels([ @headers ]);
      next;
    }

    my @j = split ",", $junction;
    my $js = parse_junction($j[0]);
    my $je = parse_junction($j[1]);
    my $known = $known_mode ? $row[2] : undef;

    if ($restrict_ref and $known eq "known") {
      die "can't find transcript header" if $transcript_index == -1;
      my %t = map {$_, 1} split /,/, $row[$transcript_index];
      next unless $t{$restrict_ref};
    }

    die if $je->{base} < $js->{base};
    die unless $js->{chrom} eq $je->{chrom};

    my %genes;

    foreach my $j ($js, $je) {
      my $base = $j->{base};
      if ($GA->find(
		    "-reference" => $j->{chrom},
		    "-start" => $base,
		    "-end" => $base
		   )) {
	my $hits = $GA->results_genes();
	foreach (@{$hits}) {
	  $genes{$_} = 1;
	}
      }
    }

    $multi_genes++ if scalar keys %genes > 1;

    next if $restrict_gene and !$genes{$restrict_gene};

    if ($preserve_mode) {
      #
      # for known junctions, preserve transcript annotations.
      # replace gene symbols as these may not be HUGO symbols for
      # all databases.
      #
      my %r;
      @r{@headers} = @row;
      $r{genes} = join ",", sort keys %genes if %genes;
      # if we don't have any symbols, leave any existing
      $rpt->end_row(\%r);
    } elsif (%genes) {
      my %r;
      $r{genes} = join ",", sort keys %genes;
      $r{junction} = $junction;
      $r{count} = $count;
      $r{known} = $known if $known_mode;
      $rpt->end_row(\%r);
    } else {
      $no_gene++;
    }
  }
  $rpt->finish();
  printf STDERR "rows skipped (no gene match): %d\n", $no_gene;
  printf STDERR "rows hitting > 1 gene: %d\n", $multi_genes;
}

sub parse_junction {
  my ($j) = @_;
  my ($chr, $pos, $strand) = split /:/, $j;
  my %r;
  $r{chrom} = $chr;
  $r{base} = $pos;
  $r{strand} = $strand;
  return \%r;
}
