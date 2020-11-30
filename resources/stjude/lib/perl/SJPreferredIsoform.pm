package SJPreferredIsoform;
# parse/retrieve SJ preferred isoforms, e.g.
# /nfs_exports/apps/gnu-apps/NextGen/nextgensupport2/NCBI/RefSeq/gene_transcript_matrix.withNM.mod

use strict;

use Configurable;
use Exporter;
use WorkingFile;

@SJPreferredIsoform::ISA = qw(Configurable Exporter);
@SJPreferredIsoform::EXPORT_OK = qw();

use MethodMaker qw(
	file
gene2preferred
nm2preferred
gene2nms
genes_ordered
gene2lines
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self) = @_;
  my $file = $self->file || die "-file";
  open(SJPI, $file) || die "can't open $file: $!";
  my %preferred;
  my %by_nm;
  my %gene2nms;
  my @genes_ordered;
  my %gene2lines;
  while (my $line = <SJPI>) {
    chomp $line;
    $line =~ s/\r$//;
    my @f = split /\t/, $line;
    die unless @f == 3;
    my ($gene, $db_id, $nm) = @f;
    my $primary;
    if ($preferred{$gene}) {
      $primary = 0;
    } else {
      $primary = 1;
      $preferred{$gene} = $nm;
      push @genes_ordered, $gene;
    }
    $by_nm{$nm} = $primary;
    push @{$gene2nms{$gene}}, $nm;
    push @{$gene2lines{$gene}}, $line;
  }
  close SJPI;
  $self->gene2preferred(\%preferred);
  $self->nm2preferred(\%by_nm);
  $self->gene2nms(\%gene2nms);
  $self->gene2lines(\%gene2lines);
  $self->genes_ordered(\@genes_ordered);
}

sub get_preferred_isoform {
  my ($self, $gene) = @_;
  return $self->gene2preferred()->{$gene};
}

sub is_preferred_nm {
  my ($self, $nm) = @_;
  return $self->nm2preferred()->{$nm};
}

sub get_nm_count {
  my ($self, $gene) = @_;
  return scalar @{$self->gene2nms->{$gene}};
}

sub get_accessions {
  my ($self, $gene) = @_;
  return $self->gene2nms->{$gene};
}

sub write_preferred_file {
  my ($self, %options) = @_;
  my $reorder = $options{"-reorder"} || {};
  my $outfile = "sjpi.tab";

  my $wf = new WorkingFile($outfile);
  my $fh = $wf->output_filehandle();
  foreach my $gene (@{$self->genes_ordered}) {
    my $lines = $self->gene2lines->{$gene} || die;
    my @postponed;
    my $promote = $reorder->{$gene};
    my @write;
    if ($promote) {
      my ($first, @rest);
      foreach my $l (@{$lines}) {
	my @f = split /\t/, $l;
	my $acc = $f[2];
	if ($acc eq $promote) {
	  $first = $l;
	} else {
	  push @rest, $l;
	}
      }
      die unless $first and @rest;
      @write = ($first, @rest);
#      die join "\n", @{$lines}, "---", @write;
    } else {
      @write = @{$lines};
    }

    foreach my $line (@write) {
      printf $fh "%s\n", $line;
    }

  }
  $wf->finish;
}

sub get_genes {
  my ($self) = @_;
  return [ sort keys %{$self->gene2nms} ];
}

sub override_preferred_isoform {
  my ($self, %options) = @_;
  my $gene = $options{"-gene"} || die;
  my $wanted_nm = $options{"-accession"} || die;
  my $all_nm = $self->gene2nms()->{$gene} || die "no accessions for $gene";
  die "can't find $wanted_nm in list for $gene" unless grep {$_ eq $wanted_nm} @{$all_nm};

  foreach my $nm (@{$all_nm}) {
    my $is_preferred = $nm eq $wanted_nm ? 1 : 0;
    $self->nm2preferred->{$nm} = $is_preferred;
  }

}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
