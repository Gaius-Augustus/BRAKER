#!/usr/bin/perl

# Author: Katharina Hoff
# E-Mail: katharina.hoff@uni-greifswald.de
# Date: May 16th 2018

# This script was developed on request of a braker user who had generated a asn
# file with ProSplign, converted that to gff-format and wanted to summarize
# multiplicity information over that gff-file for intron hints. This output
# is compatible with braker.pl, but not directly with gmes_petap.pl (GeneMark-EX).

# Input format example:
# seq1	ProSplign	exon	12839	13668	.	-	.	grp=prot1;pri=4;src=P
# seq1	ProSplign	intron	12752	12835	.	-	.	grp=prot1;pri=4;src=P
# seq1	ProSplign	intron	12752	12835	.	-	.	grp=prot1;pri=4;src=P

# Output format example:
# seq1	ProSplign	intron	12752	12835	.	-	.	mult=1;pri=4;src=P
# seq1	ProSplign	intron	12752	12835	.	-	.	mult=1;pri=4;src=P

# usage: cat hintsfile | summarize_hints.pl > summarized_hintsfile

my %hints;

while(<STDIN>){
	if(m/\tintron\t/){
		my @f = split(/\t/);
		my $hint = "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\t";
		if(not(defined($hints{$hint}))) {
			$hints{$hint} = 1;
		}else{
			$hints{$hint} += 1;
		}
	}
}

foreach( keys %hints) {
	print $_."mult=$hints{$_};src=P;pri=4;\n"
}