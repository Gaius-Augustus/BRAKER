#!/usr/bin/perl
# ==============================================================
# Katharina J. Hoff
# This script takes a hints file generated from the GATECH
# protein mapping pipeline (contains CDSpart, start and stop
# hints) and outputs a list of single exon genes that have
# start and stop
#
# last changes: February 13th 2018
# ==============================================================

use strict;
use warnings;

use Getopt::Long qw( GetOptions );
use FindBin qw( $RealBin );
use File::Spec;
use Cwd qw( abs_path cwd );
use Data::Dumper;
use YAML;

# ------------------------------------------------
my $v = 0;
my $debug = 0;

my $cfg;
my $log;

my $bin = $RealBin;
my $work_dir = cwd;
# ------------------------------------------------
my $in_gff_file  = '';
my $out_gff_file = '';
# ------------------------------------------------

Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();

# ------------------------------------------------
sub CheckBeforeRun
{
	print "check before run\n" if $debug;
		
	$bin      = ResolvePath( $bin );
	$work_dir = ResolvePath( $work_dir );
	
	$in_gff_file = ResolvePath( $in_gff_file );
	
	if( !$in_gff_file )  { print "error, required file name is missing $0:  option --in_gff\n"; exit 1; }
	if( !$out_gff_file ) { print "error, required file name is missing $0:  option --out_gff\n"; exit 1; }
};
# ------------------------------------------------
sub ResolvePath
{
	my( $name, $path ) = @_;
	return '' if !$name;
	$name = File::Spec->catfile( $path, $name ) if ( defined $path and $path );
	if( ! -e $name ) { print "error, file not found $0: $name\n"; exit 1; }
	return abs_path( $name );
};
# ------------------------------------------------
sub ParseCMD
{
	print "parse cmd\n" if $debug;
	
	my $cmd = $0;
	foreach my $str (@ARGV) { $cmd .= ( ' '. $str ); }
	
	my $opt_results = GetOptions
	(
		'in_gff=s'  => \$in_gff_file,
		'out_gff=s' => \$out_gff_file,
		'verbose' => \$v,
		'debug'   => \$debug
	);

	if( !$opt_results ) { print "error on command line: $0\n"; exit 1; }
	if( @ARGV > 0 ) { print "error, unexpected argument found on command line: $0 @ARGV\n"; exit 1; }
	$v = 1 if $debug;

	# save informaton for debug
	$cfg->{'d'}->{'in_gff_file'}  = $in_gff_file;
	$cfg->{'d'}->{'out_gff_file'} = $out_gff_file;
	$cfg->{'d'}->{'v'}     = $v;
	$cfg->{'d'}->{'debug'} = $debug;
	$cfg->{'d'}->{'cmd'}   = $cmd;
	
#	print Dumper($cfg) if $debug;
};
# ------------------------------------------------
sub Usage
{
	print qq(# -------------------
Usage:  $0

Required options:
  --in_gff   [name] name of file with ProSplign hints
  --out_gff  [name] output
 
Developer options:
  --verbose
  --debug
# -------------------
);
	exit 1;
}
# ================== END sub =====================