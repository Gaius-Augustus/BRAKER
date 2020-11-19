package Test::BRAKER;

use strict;
use warnings;
use Exporter qw(import);
use IO::File   ();
use Test::More ();

our @EXPORT_OK = qw(fasta_file_ok);

sub fasta_file_ok {
  my ($desc, @files) = (pop(@_), shift, shift);
  my @recordsets;

  for my $i (0 .. 1) {
    my $fh = IO::File->new($files[$i]);
    local $/ = ">";
    my $first = <$fh>;
    while (<$fh>) {
      chomp;
      my ($header, $record) = split /\n/, $_, 2;
      $record =~ s/\s//gs;
      $recordsets[$i]{$header} = $record;
    }
  }

  my ($results);
  for my $id (sort keys %{$recordsets[0]}) {
    $results += _test('is', $recordsets[0]{$id}, $recordsets[1]{$id}, "sequence for $id match");
  }

  return _test('is', $results, scalar(keys %{$recordsets[1]}), $desc);
}

sub _test {
  my ($method_name, @args) = @_;
  local $Test::Builder::Level = $Test::Builder::Level + 3;
  return Test::More->can($method_name)->(@args);
}

1;

__END__

=encoding utf8

=head1 NAME

Test::BRAKER - Helper module for testing BRAKER

=head1 SYNOPSIS

  use Test::More;
  use lib 't/lib';
  use Test::BRAKER qw(fasta_file_ok);

  fasta_file_ok($fasta_file, $fixture_fasta, 'generated fasta file congruent');

=head1 DESCRIPTION

Some L<Test::More> style functions to aid testing parts of L<BRAKER|https://github.com/Gaius-Augustus/BRAKER>.

=cut
