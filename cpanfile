# -*- mode: perl -*-

requires 'Exporter';          # core since 5
requires 'File::Basename';    # core since 5
requires 'File::Spec';        # core since 5.00405
requires 'Hash::Merge';
requires 'List::Util';        # core since v5.7.3
requires 'MCE::Mutex';
requires 'Math::Utils';
requires 'Module::Load::Conditional';
requires 'Parallel::ForkManager';
requires 'POSIX';             # core since 5
requires 'Scalar::Util::Numeric';
requires 'YAML';

test_requires 'Cwd';           # core since 5
test_requires 'File::Glob';    # core since v5.6.0
test_requires 'File::Temp';    # core since v5.6.1
test_requires 'IO::File';      # core since 5.00307
test_requires 'Test::More';    # core since v5.6.2

feature prothint => sub {
  requires 'threads';          # core since v5.7.3
};
