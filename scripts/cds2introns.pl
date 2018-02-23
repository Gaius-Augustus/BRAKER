#!/usr/bin/perl

# Author: Katharina J. Hoff, Feb 21st 2018

# This script extracts introns in gtf format from a gtf CDS file.

if(scalar(@ARGV) != 1){
	die("exons2introns.pl exons.gtf > introns.gtf\n\nThis script extracts introns in gtf format from a gtf exon file.\n\nIt is assumed that the exon file was sorted with\n\nsort -k4,4 -g | sort -k9,9\n\n")
}

my %genes = ();
my @t;

# read exon gtf file
open(EXONS, "<", $ARGV[0]) or die "Could not open exon file $ARGV[0]!\n";
while(<EXONS>){
	if(m/\tCDS\t/){
		@t = split(/\t/);
		push(@{$genes{$t[8]}}, $_);
	}
}
close(EXONS) or die "Could not close exon file $ARGV[0]!\n";

my $ele;
for $transcript (keys %genes){
	$ele = 1;
	foreach(@{$genes{$transcript}}){
		@t = split(/\t/);
		if($ele > 1){
			print "$t[0]\t$t[1]\tintron\t$intronStart\t".($t[3]-1)."\t0\t$t[6]\t$t[7]\t$t[8]";
		}
		$intronStart = $t[4]+1;
		$ele = $ele + 1;
	}
}

