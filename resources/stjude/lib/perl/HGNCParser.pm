package HGNCParser;
# parser/lookup for HGNC (HUGO) gene symbols & synonyoms

use strict;
use Configurable;
use Exporter;

use DelimitedFile;
use MiscUtils qw(dump_die);

@HGNCParser::ISA = qw(Configurable Exporter);
@HGNCParser::EXPORT_OK = qw();

use MethodMaker qw(
	file

	ids

	index_approved
	index_refseq
	index_synonym
	index_prev
	index_all_sym
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
  my ($self, %options) = @_;
  my $fn = $self->file || die;
  my $df = new DelimitedFile(
			     "-file" => $fn,
			     "-headers" => 1,
			     );

  my %ids;

  my %idx_approved;
  my %idx_nm;
  my %idx_synonym;
  my %idx_prev;
  my %idx_all_sym;

  my %allsym2id;
  # all symbols

  while (my $row = $df->get_hash()) {
    my $id = $row->{"HGNC ID"} || die;
    die "duplicate $id" if $ids{$id};
    $ids{$id} = $row;
    # store rows by unique ID, indexes point to IDs

    my $approved_sym = $row->{"Approved Symbol"} || die;
    $idx_approved{$approved_sym}{$id} = 1;
    $idx_approved{uc($approved_sym)}{$id} = 1;
    # use same structure as other indexes

    my @all_syms = $approved_sym;

    if (my $refs = $row->{"RefSeq IDs"}) {
      my @rs = split /,\s*/, $refs;
      foreach my $rs (@rs) {
	$idx_nm{$rs}{$id} = 1;
      }
      $row->{refseqs_hash} = { map {$_, 1} @rs };
    }

    foreach my $ref (
		     [ "Previous Symbols", \%idx_prev ],
		     [ "Synonyms", \%idx_synonym ],
		    ) {
      my ($field, $hash) = @{$ref};
      die unless exists $row->{$field};
      if (my $things = $row->{$field}) {
	my @things = split /,\s*/, $things;
	push @all_syms, @things;
	foreach my $thing (@things) {
	  $hash->{uc($thing)}{$id} = 1;
	}
      }
    }
#    printf STDERR "all for %s => %s\n", $approved_sym, join ",", @all_syms;

    foreach my $sym (@all_syms) {
      $idx_all_sym{uc($sym)}{$id} = 1;
    }
  }

  $self->ids(\%ids);
  $self->index_approved(\%idx_approved);
  $self->index_refseq(\%idx_nm);
  $self->index_prev(\%idx_prev);
  $self->index_synonym(\%idx_synonym);
  $self->index_all_sym(\%idx_all_sym);
}

sub find {
  my ($self, %options) = @_;
  my $symbol = $options{"-symbol"};
  my @ids;
  if ($symbol) {
    $symbol = uc($symbol);
    my $hit;
    if ($options{"-previous"}) {
      $hit = $self->find_by_index($self->index_prev(), $symbol);
    } elsif ($options{"-approved"}) {
      $hit = $self->find_by_index($self->index_approved(), $symbol);
    } elsif ($options{"-synonym"}) {
      $hit = $self->find_by_index($self->index_synonym(), $symbol);
    } else {
      # any symbol
      $hit = $self->find_by_index($self->index_all_sym(), $symbol);
    }
    die if $hit and not(ref $hit);
    @ids = keys %{$hit} if $hit;
  } else {
    # other types?
    die "-symbol";
  }

  my @results;
  if (@ids) {
    @results = map {$self->ids->{$_} || die "no ID for $_"} @ids;
  }

  return @results ? \@results : undef;
}

sub is_approved {
  my ($self, $sym) = @_;
  my $status = $self->index_approved->{$sym};
  printf STDERR "status for %s is %s\n", $sym, $status;
  return $status;
}

sub find_by_index {
  my ($self, $index, $sym) = @_;
  return $index->{uc($sym)};
  # symbols are always uppercased in index
}

sub prune_synonyms {
  my ($self) = @_;

  my $idx_synonym = $self->index_synonym();
  my $idx_approved = $self->index_approved();
  foreach my $sym (sort keys %{$idx_synonym}) {
    my $bad;
    if ($idx_approved->{$sym}) {
      # if a synonym is itself an approved symbol, remove link to other gene.
      # e.g. FAH is both an approved symbol as well as a synonym for FANCA.
      $bad = 1;
    }

    if (keys %{$idx_synonym->{$sym}} > 1) {
      # e.g. 2F1 is ambiguous: maps to both SLC25A5 and KLRG1
#      printf STDERR "ambiguous synonym %s => %s\n", $sym, join ",", sort keys %{$idx_synonym->{$sym}};
      $bad = 1;
    }

    delete $idx_synonym->{sym} if $bad;
  }

}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
