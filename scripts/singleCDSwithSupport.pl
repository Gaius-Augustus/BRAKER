#!/usr/bin/perl
# ==============================================================
# Katharina J. Hoff and Alexandre Lomsadze
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
my $v     = 0;
my $debug = 0;

my $cfg;
my $log;

my $bin      = $RealBin;
my $work_dir = cwd;

# ------------------------------------------------
my $in_gff_file              = '';
my $out_gff_file             = '';
my $stopCodonExcludedFromCDS = 0;

# ------------------------------------------------

my %h;
my %s;

Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();
ReadGff();
SummarizeHints();
PrintCDS();

# ------------------------------------------------
sub PrintCDS {
    open( my $OUT, ">", $out_gff_file )
        or die("$!, error on open file $out_gff_file");
    foreach my $txid ( keys %s ) {
        print $OUT $s{$txid}{'seq'} . "\t"
            . $s{$txid}{'src'}
            . "\tCDS\t"
            . $s{$txid}{'start'} . "\t";
        if ( $stopCodonExcludedFromCDS == 0 ) {
            print $OUT $s{$txid}{'stop'};
        }
        else {
            print $OUT ( $s{$txid}{'stop'} - 3 );
        }
        print $OUT "\t"
            . $s{$txid}{'score'} . "\t"
            . $s{$txid}{'strand'} . "\t"
            . $s{$txid}{'frame'} . "\t";
        if ( $s{$txid}{'mult'} > 1 ) {
            print $OUT "src=P;pri=4;mult=" . $s{$txid}{'mult'} . ";\n";
        }
        else {
            print $OUT $s{$txid}{'grp'} . "\n";
        }
    }
    close($OUT) or die("$!, error on close file $out_gff_file");
}

# ------------------------------------------------
sub SummarizeHints {
    foreach my $txid ( keys %h ) {
        if (   defined( $h{$txid}{'start'} )
            && defined( $h{$txid}{'stop'} )
            && !defined( $h{$txid}{'introns'} )
            && $h{$txid}{'cds'} == 1 )
        {
            my $cdsid
                = $h{$txid}{'seq'} . "_"
                . $h{$txid}{'start'} . "_"
                . $h{$txid}{'stop'} . "_"
                . $h{$txid}{'strand'};
            if ( !defined( $s{$cdsid} ) ) {
                $s{$cdsid}{'seq'}    = $h{$txid}{'seq'};
                $s{$cdsid}{'src'}    = $h{$txid}{'src'};
                $s{$cdsid}{'start'}  = $h{$txid}{'start'};
                $s{$cdsid}{'stop'}   = $h{$txid}{'stop'};
                $s{$cdsid}{'score'}  = $h{$txid}{'score'};
                $s{$cdsid}{'strand'} = $h{$txid}{'strand'};
                $s{$cdsid}{'frame'}  = $h{$txid}{'frame'};
                $s{$cdsid}{'grp'}    = $h{$txid}{'grp'};
                $s{$cdsid}{'mult'}   = 1;

            }
            else {
                $s{$cdsid}{'mult'}++;
            }
        }
    }
}

# ------------------------------------------------
sub ReadGff {
    open( my $IN, "<", $in_gff_file )
        or die("$!, error on open file $in_gff_file");
    while (<$IN>) {
        $_ =~ m/grp=(.*);/;
        my $txid = $1;
        if ( $_ =~ m/\tintron\t/ ) {
            if ( !defined( $h{$txid}{'introns'} ) ) {
                $h{$txid}{'introns'} = 1;
            }
            else {
                $h{$txid}{'introns'}++;
            }

        }
        elsif ( $_ =~ m/\tCDS/ ) {
            if ( !defined( $h{$txid}{'cds'} ) ) {
                $h{$txid}{'cds'} = 1;
            }
            else {
                $h{$txid}{'cds'}++;
            }
        }
        elsif ( $_
            =~ m/^(.*)\t(.*)\tstart\t(\d+)\t(\d+)\t(.*)\t(\+|-)\t(.*)\t(.*)/ )
        {
            $h{$txid}{'seq'}    = $1;
            $h{$txid}{'src'}    = $2;
            $h{$txid}{'score'}  = $5;
            $h{$txid}{'frame'}  = $7;
            $h{$txid}{'strand'} = $6;
            if ( $h{$txid}{'strand'} eq '+' ) {
                $h{$txid}{'start'} = $3;
            }
            else {
                $h{$txid}{'stop'} = $4;
            }
            $h{$txid}{'grp'} = $8;
        }
        elsif ( $_ =~ m/^(.*)\t(.*)\tstop\t(\d+)\t(\d+)\t(.*)\t(\+|-)\t/ ) {
            if ( $5 eq '+' ) {
                $h{$txid}{'stop'} = $4;
            }
            else {
                $h{$txid}{'stop'} = $3;
            }

        }
    }
    close($IN) or die("$!, error on close file $in_gff_file");
}

# ------------------------------------------------
sub CheckBeforeRun {
    print "check before run\n" if $debug;

    $bin      = ResolvePath($bin);
    $work_dir = ResolvePath($work_dir);

    $in_gff_file = ResolvePath($in_gff_file);

    if ( !$in_gff_file ) {
        print "error, required file name is missing $0:  option --in_gff\n";
        exit 1;
    }
    if ( !$out_gff_file ) {
        print "error, required file name is missing $0:  option --out_gff\n";
        exit 1;
    }
}

# ------------------------------------------------
sub ResolvePath {
    my ( $name, $path ) = @_;
    return '' if !$name;
    $name = File::Spec->catfile( $path, $name )
        if ( defined $path and $path );
    if ( !-e $name ) { print "error, file not found $0: $name\n"; exit 1; }
    return abs_path($name);
}

# ------------------------------------------------
sub ParseCMD {
    print "parse cmd\n" if $debug;

    my $cmd = $0;
    foreach my $str (@ARGV) { $cmd .= ( ' ' . $str ); }

    my $opt_results = GetOptions(
        'in_gff=s'                   => \$in_gff_file,
        'out_gff=s'                  => \$out_gff_file,
        'stopCodonExcludedFromCDS:i' => \$stopCodonExcludedFromCDS,
        'verbose'                    => \$v,
        'debug'                      => \$debug
    );

    if ( !$opt_results ) { print "error on command line: $0\n"; exit 1; }
    if ( @ARGV > 0 ) {
        print "error, unexpected argument found on command line: $0 @ARGV\n";
        exit 1;
    }
    $v = 1 if $debug;

    # save informaton for debug
    $cfg->{'d'}->{'in_gff_file'}              = $in_gff_file;
    $cfg->{'d'}->{'out_gff_file'}             = $out_gff_file;
    $cfg->{'d'}->{'v'}                        = $v;
    $cfg->{'d'}->{'debug'}                    = $debug;
    $cfg->{'d'}->{'cmd'}                      = $cmd;
    $cfg->{'d'}->{'stopCodonExcludedFromCDS'} = $stopCodonExcludedFromCDS;

    #	print Dumper($cfg) if $debug;
}

# ------------------------------------------------
sub Usage {
    print qq(# -------------------
Usage:  $0

Required options:
  --in_gff                     [name] name of file with ProSplign hints
  --out_gff                    [name] output

Developer options:
  --stopCodonExcludedFromCDS=1 default 0
  --verbose
  --debug
# -------------------
);
    exit 1;
}

# ================== END sub =====================
