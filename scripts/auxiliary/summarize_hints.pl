#!/usr/bin/perl

# Author: Katharina Hoff
# E-Mail: katharina.hoff@uni-greifswald.de
# Date: May 16th 2018

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