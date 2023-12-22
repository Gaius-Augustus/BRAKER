# -*- mode: perl; -*-
use strict;
use warnings;
use Test::More;
use Cwd ();
use File::Basename qw(basename);
use File::Glob ':bsd_glob';
use File::Temp qw(tempdir);

our $exit_called = 0;

BEGIN {
  # this override needs to be in a begin block
  # see ikegami's advice https://stackoverflow.com/a/25381153
  *CORE::GLOBAL::exit = sub (;$) { $exit_called++ };
}

use lib 'scripts';
use helpMod_braker;

subtest 'Directory removed by clean_abort' => sub {
  helpMod_braker->import('clean_abort');

  my $dir = tempdir('util.XXXXX', TMPDIR => 1);
  open my $fh, '>', File::Spec->catfile($dir, 'test.cfg');
  close $fh;

  my $buffer = '';
  open my $handle, '>', \$buffer;
  local *STDERR = $handle;

  local $exit_called = 0;
  clean_abort($dir, 0, 'simple message');

  is $exit_called, 1,                'exit called';
  is $buffer,      'simple message', 'message reported to stderr';
  ok !-d $dir, 'removed';
};

subtest 'Directory is not removed by clean_abort' => sub {
  helpMod_braker->import('clean_abort');

  my $dir = tempdir('util.XXXXX', TMPDIR => 1, CLEANUP => 1);
  open my $fh, '>', File::Spec->catfile($dir, 'test.cfg');
  close $fh;

  my $buffer = '';
  open my $handle, '>', \$buffer;
  local *STDERR = $handle;

  local $exit_called = 0;
  clean_abort($dir, 1, 'simple message');

  is $exit_called, 1,                'exit called';
  is $buffer,      'simple message', 'message reported to stderr';
  is -d $dir, 1, 'directory still exists, has not been removed';
};


done_testing;
