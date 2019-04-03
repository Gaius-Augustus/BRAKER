#!/usr/bin/env perl

####################################################################################################
#                                                                                                  #
# filter_augustus_gff.pl                                                                           #
#                                                                                                  #
# Author: Katharina Hoff                                                                           # 
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
####################################################################################################

use strict;
use warnings;
use Getopt::Long;

my $usage = <<'ENDUSAGE';

filter_augustus_gff.pl     Filter native AUGUSTUS output for genes that have evidential support

SYNOPSIS

filter_augustus_gff.pl --in=augustus.gff --out=filtered.gtf --src=(E|P)

INPUT FILE OPTIONS

--in=augustus.gff                   AUGUSTUS output file in GFF format
--out=filtered.gff                  If specified, write filtered gene predictions
                                    in GTF format to given file. If not specified,
                                    only the number of supported genes is printed
                                    to STDOUT.
--src=(E|P)                         Source tag of evidence to filter for, e.g. P 
                                    for protein or E for ESTs/RNA-Seq introns,
                                    can in principle be any source, but only one 
                                    at a time. Support threshold is one hint per 
                                    transcript.
--version                           Print version number of filter_augustus_gff.pl
--help                              Print this help message

DESCRIPTION

Examples:

filter_augustus_gff.pl --in=augustus.gff --out=filtered.gtf --src=P 

filter_augustus_gff.pl --in=augustus.gff --src=P

ENDUSAGE

my $version = "1.0";
my $vers;
my $help;
my $in_file;
my $out_file;
my $src;

if ( @ARGV == 0 ) {
    print "$usage\n";
    exit(0);
}

GetOptions(
    'in=s'     => \$in_file,
    'out=s'    => \$out_file,
    'src=s'    => \$src,
    'version!' => \$vers,
    'help!'    => \$help);

my @gtf_lines;
my $tx_id;
my $supported_by;
my $thereof_src;
my $start_screening = 0;
my %keep;

if($help){
	print $usage;
	exit(0);
}

if($vers){
	print "version ".$version."\n";
	exit(0);
}

if(not(defined($in_file))){
	print "ERROR: in file " . __FILE__ ." at line "
	      . __LINE__ ."\nNo input file provided!\n";
    print $usage;
    exit(1);
}

if(not(defined($src))){
	print "ERROR: in file " . __FILE__ ." at line "
          . __LINE__ ."\nNo src for filtering provided!\n";
    print $usage;
    exit(1);
}


open(INFILE, "<", $in_file) or die ("ERROR: in file " . __FILE__ ." at line "
                                     . __LINE__ ."\nFailed to open file $in_file for reading!\n");

while(<INFILE>){
	if(m/transcript_id \"([^;]+)\";/){
		$tx_id = $1;
	}elsif(m/hint groups fully obeyed:\s+(\d+)/){
		$supported_by = $1;
		$start_screening = 1;
		$thereof_src = 0;
	}elsif(m/incompatible hint groups:/){
		$start_screening = 0;
	}elsif(m/$src:\s+(\d+)/ && $start_screening==1){
		$thereof_src = $1;
		if($thereof_src > 0) {
			$keep{$tx_id} = $thereof_src;
		}
	}
	if(m/transcript_id \"/){
		push(@gtf_lines, $_);
	}
}

close(INFILE) or die ("ERROR: in file " . __FILE__ ." at line "
                                     . __LINE__ ."\nFailed to close $in_file!\n");

if(defined($out_file)){
	open(OUTFILE, ">", $out_file) or die ("ERROR: in file " . __FILE__ ." at line "
                                          . __LINE__ ."\nFailed to open file $out_file for writing!\n");
	foreach(@gtf_lines){
		if(m/transcript_id \"([^;]+)\";/){
			if(defined($keep{$1})){
				print OUTFILE $_;
			}
		}
}
	close(OUTFILE) or die ("ERROR: in file " . __FILE__ ." at line "
                                      . __LINE__ ."\nFailed to close $out_file!\n");
}else{
	my $n_tx = keys %keep;
	print "$n_tx\n";
}