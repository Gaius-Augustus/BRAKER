# -*- mode: perl; -*-
use strict;
use warnings;
use Test::More;
use Cwd         ();
use File::Basename qw(basename);
use File::Glob ':bsd_glob';
use File::Temp qw(tempdir);
use IO::File ();

use lib 't/lib';
use Test::BRAKER qw(fasta_file_ok);
no lib 't/lib';

BEGIN {
  use lib 'scripts';
  use_ok 'helpMod_braker';
}

subtest 'checkFile' => sub {
  helpMod_braker->import('checkFile');

  # absolute
  my $file = eval { checkFile(Cwd::abs_path('scripts/helpMod_braker.pm'), 'perl module file'); };
  is $file, Cwd::abs_path('scripts/helpMod_braker.pm'), 'absolute path found';
  is $@, '', 'exists, no error';

  # relative
  $file = eval { checkFile('scripts/helpMod_braker.pm', 'perl module file'); };
  is $file, Cwd::abs_path('scripts/helpMod_braker.pm'), 'path now absolute';
  is $@, '', 'exists so no error';

  # incorrect call
  eval { checkFile('', 'unknown', 'usage:') };
  like $@, qr{^ERROR: in file .*helpMod_braker\.pm at line [0-9]+}, 'match message';
  like $@, qr{^missing unknown file}m,                       'missing input';

  # tilde expansion
  my $dir    = tempdir('util.XXXXX', TMPDIR => 1, CLEANUP => 1);
  my $subdir = File::Spec->catfile($dir, 'test-glob');
  mkdir $subdir or fail();
  my $path = File::Spec->catfile($subdir, 'script.pl');
  open my $fh, '>', $path or fail();
  $file = eval {
    local $ENV{HOME} = $dir;
    my $tilde_path = File::Spec->catfile('~', basename($subdir), 'script.pl');
    checkFile($tilde_path, 'perl script');
  };
  is $@, '', 'no error';
  is $file, $path, 'correct path expansion';

  # absolute path that does not exist/readable
  eval { checkFile('/root/root/root/root/does/not/exist', 'binary', 'test-suite-eval'); };
  like $@, qr{^ERROR: in file [^\s]+ at line [0-9]+},                        'match message';
  like $@, qr{^ binary file /root/root/root/root/does/not/exist not found}m, 'not found';
};

subtest 'find' => sub {
  helpMod_braker->import('find');
  my $result = find('braker.pl', 'scripts', 'scripts', 'scripts');
  ok $result;
  is basename($result), 'braker.pl', 'correct script';
};

subtest 'formatDetector' => sub {
  helpMod_braker->import('formatDetector');
  my $fmt = formatDetector('example/genome.fa');
  {
    local $TODO = 'genomes are nucleotides';
    is $fmt, 'fasta-dna', 'correct';
  }

  $fmt = formatDetector('example/proteins.fa');
  is $fmt, 'fasta-prot', 'correct';

  $fmt = formatDetector('example/RNAseq.hints');
  is $fmt, 'gff', 'hints file is gff';

  $fmt = formatDetector('example/results/test1/augustus.hints.gtf');
  is $fmt, 'gff', 'hints file is gff';

  $fmt = formatDetector('example/results/test1/genemark.gtf');
  is $fmt, 'gff', 'genemark file appears as gff';

  $fmt = formatDetector('t/data/sequence.gb');
  is $fmt, 'gb', 'fixture is in genbank format';

  # missing file detection
  $fmt = eval { formatDetector('example/does-not-exist.gff') };
  like $@, qr{^Could not open example/does-not-exist\.gff}m, 'error message match';
  is $fmt, undef,                                            'no format returned';

  # damaged fasta format detection
  my $dir        = tempdir('fmt-detect.XXXXX', TMPDIR => 1, CLEANUP => 1);
  my $damaged_fa = File::Spec->catfile($dir, 'damaged.fa');
  my $buffer     = '';
  {
    my $fh = IO::File->new($damaged_fa, '>');
    print $fh join "\n", ('>') x 100;
    close $fh;    # close to flush to disk

    open my $handle, '>', \$buffer;
    local *STDERR = $handle;
    $fmt = formatDetector($damaged_fa);
  }
  like $buffer, qr{^ERROR: in file},                                     'message damaged fasta format';
  like $buffer, qr{^$damaged_fa appears to be in corrupt FASTA format}m, 'message to user correct';
  is $fmt,      '',                                                      'no format returned';

};

subtest 'gtf2fasta' => sub {
  helpMod_braker->import('gtf2fasta');
  my $gtf_file = 'example/results/test1/augustus.hints.gtf';
  my $genome   = 'example/genome.fa';
  my $dir      = tempdir('fasta.XXXXX', TMPDIR => 1, CLEANUP => 1);
  my $fasta    = File::Spec->catfile($dir, 'coding.aa');
  gtf2fasta($genome, $gtf_file, $fasta, 1);
  ok -f $fasta, 'fasta file exists';
  fasta_file_ok($fasta, 'example/results/test1/augustus.hints.aa', 'all records match');
};

subtest 'reverse_complement' => sub {
  is helpMod_braker::reverse_complement('AGGGGGGGCCCTTA'), 'TAAGGGCCCCCCCT', 'correct';
  is helpMod_braker::reverse_complement('ATAT'),           'ATAT',           'correct';
  is helpMod_braker::reverse_complement('ATATa'),          'tATAT',          'correct';
  is helpMod_braker::reverse_complement('atgc'),           'gcat',           'correct';

  # incorrect input
  is helpMod_braker::reverse_complement(''),    '',     'correct';
  is helpMod_braker::reverse_complement('ZXS'), 'SXZ',  'odd, but correct';
  is helpMod_braker::reverse_complement(1000),  '0001', 'odd, but correct';
};

subtest 'tildeConvert' => sub {
  helpMod_braker->import('tildeConvert');
  my $dir = tempdir('util.XXXXX', TMPDIR => 1, CLEANUP => 1);
  mkdir File::Spec->catfile($dir, 'test-glob');
  for my $path_to_glob (qw(~/ ~ ~/test-glob)) {
    local $ENV{HOME} = $dir;
    is tildeConvert($path_to_glob), bsd_glob($path_to_glob, GLOB_TILDE | GLOB_ERR), 'correct';
  }
};

subtest 'tildeConvert edge case' => sub {
  my $homedir = bsd_glob('~sudo', GLOB_TILDE | GLOB_ERR);
  is $homedir, undef, 'no users are called sudo - sanity';
TODO: {
    local $TODO = '~sudo the the home of sudo user';
    is tildeConvert('~sudo/augustus'), bsd_glob('~sudo/augustus', GLOB_TILDE | GLOB_ERR), 'future';
  }
};

subtest 'Translate with dna2aa' => sub {

  my @codons = qw(
    ATGGGA
    ATGTAA
    ATGTAG
    ATGTGA
    ATGCTGTAA
  );

  # initial reference (1)
  my $code = 1;
  is helpMod_braker::dna2aa($codons[0], $code), 'MG',  'correct GGA';
  is helpMod_braker::dna2aa($codons[1], $code), 'M*',  'correct TAA';
  is helpMod_braker::dna2aa($codons[2], $code), 'M*',  'correct TAG';
  is helpMod_braker::dna2aa($codons[3], $code), 'M*',  'correct TGA';
  is helpMod_braker::dna2aa($codons[4], $code), 'ML*', 'correct CTG TAA';

  $code = 6;
  is helpMod_braker::dna2aa($codons[0], $code), 'MG',  'correct GGA';
  is helpMod_braker::dna2aa($codons[1], $code), 'MQ',  'correct TAA';
  is helpMod_braker::dna2aa($codons[2], $code), 'MQ',  'correct TAG';
  is helpMod_braker::dna2aa($codons[3], $code), 'M*',  'correct TGA';
  is helpMod_braker::dna2aa($codons[4], $code), 'MLQ', 'correct CTG TAA';

  $code = 10;
  is helpMod_braker::dna2aa($codons[0], $code), 'MG',  'correct GGA';
  is helpMod_braker::dna2aa($codons[1], $code), 'M*',  'correct TAA';
  is helpMod_braker::dna2aa($codons[2], $code), 'M*',  'correct TAG';
  is helpMod_braker::dna2aa($codons[3], $code), 'MC',  'correct TGA';
  is helpMod_braker::dna2aa($codons[4], $code), 'ML*', 'correct CTG TAA';

  $code = 12;
  is helpMod_braker::dna2aa($codons[0], $code), 'MG',  'correct GGA';
  is helpMod_braker::dna2aa($codons[1], $code), 'M*',  'correct TAA';
  is helpMod_braker::dna2aa($codons[2], $code), 'M*',  'correct TAG';
  is helpMod_braker::dna2aa($codons[3], $code), 'M*',  'correct TGA';
  is helpMod_braker::dna2aa($codons[4], $code), 'MS*', 'correct CTG TAA';

  $code = 25;
  is helpMod_braker::dna2aa($codons[0], $code), 'MG',  'correct GGA';
  is helpMod_braker::dna2aa($codons[1], $code), 'M*',  'correct TAA';
  is helpMod_braker::dna2aa($codons[2], $code), 'M*',  'correct TAG';
  is helpMod_braker::dna2aa($codons[3], $code), 'MG',  'correct TGA';
  is helpMod_braker::dna2aa($codons[4], $code), 'ML*', 'correct CTG TAA';

  $code = 26;
  is helpMod_braker::dna2aa($codons[0], $code), 'MG',  'correct GGA';
  is helpMod_braker::dna2aa($codons[1], $code), 'M*',  'correct TAA';
  is helpMod_braker::dna2aa($codons[2], $code), 'M*',  'correct TAG';
  is helpMod_braker::dna2aa($codons[3], $code), 'M*',  'correct TGA';
  is helpMod_braker::dna2aa($codons[4], $code), 'MA*', 'correct CTG TAA';

  $code = 27;
  is helpMod_braker::dna2aa($codons[0], $code), 'MG',  'correct GGA';
  is helpMod_braker::dna2aa($codons[1], $code), 'MQ',  'correct TAA';
  is helpMod_braker::dna2aa($codons[2], $code), 'MQ',  'correct TAG';
  is helpMod_braker::dna2aa($codons[3], $code), 'M*',  'correct TGA';
  is helpMod_braker::dna2aa($codons[4], $code), 'MLQ', 'correct CTG TAA';

  $code = 28;
  is helpMod_braker::dna2aa($codons[0], $code), 'MG',  'correct GGA';
  is helpMod_braker::dna2aa($codons[1], $code), 'M*',  'correct TAA';
  is helpMod_braker::dna2aa($codons[2], $code), 'M*',  'correct TAG';
  is helpMod_braker::dna2aa($codons[3], $code), 'M*',  'correct TGA';
  is helpMod_braker::dna2aa($codons[4], $code), 'ML*', 'correct CTG TAA';

  $code = 29;
  is helpMod_braker::dna2aa($codons[0], $code), 'MG',  'correct GGA';
  is helpMod_braker::dna2aa($codons[1], $code), 'MY',  'correct TAA';
  is helpMod_braker::dna2aa($codons[2], $code), 'MY',  'correct TAG';
  is helpMod_braker::dna2aa($codons[3], $code), 'M*',  'correct TGA';
  is helpMod_braker::dna2aa($codons[4], $code), 'MLY', 'correct CTG TAA';

  $code = 30;
  is helpMod_braker::dna2aa($codons[0], $code), 'MG',  'correct GGA';
  is helpMod_braker::dna2aa($codons[1], $code), 'ME',  'correct TAA';
  is helpMod_braker::dna2aa($codons[2], $code), 'ME',  'correct TAG';
  is helpMod_braker::dna2aa($codons[3], $code), 'M*',  'correct TGA';
  is helpMod_braker::dna2aa($codons[4], $code), 'MLE', 'correct CTG TAA';

  $code = 31;
  is helpMod_braker::dna2aa($codons[0], $code), 'MG',  'correct GGA';
  is helpMod_braker::dna2aa($codons[1], $code), 'M*',  'correct TAA';
  is helpMod_braker::dna2aa($codons[2], $code), 'M*',  'correct TAG';
  is helpMod_braker::dna2aa($codons[3], $code), 'MW',  'correct TGA';
  is helpMod_braker::dna2aa($codons[4], $code), 'ML*', 'correct CTG TAA';

  # dangling bases / incomplete codon
  is helpMod_braker::dna2aa('ATGGAGGAGGGTG', 1), 'MEEG', 'correct';
};


done_testing;
