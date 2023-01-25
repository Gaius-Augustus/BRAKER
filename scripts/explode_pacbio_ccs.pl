#!/usr/bin/env perl

# Blow up fastq file with assembled PacBio subreads,
# i.e. copy reads as many times as they have coverage
# according to header of read in original fastq
# This is a temporary solution, only.

# Katharina J. Hoff
# Nov 5 2021
# Contact: katharina.hoff@uni-greifswald.de

my $line = 0;
my $entry1 = "";
my $entry2 = "";
my $cov = 0;
while(<STDIN>){
	if($line == 0){
		$_ =~ m/full_length_coverage=(\d+);/;
		$cov = $1;
		if($cov == 0){
			print STDERR "Something went wrong with parsing coverage from header!\n";
			exit(1);
		}
		$entry1 = "";
		$entry2 = "";
		$entry1 .= $_;
		$line++;
	}elsif(($line == 1) or ($line == 2)){
		$entry2 .= $_;
		$line++;
	}else{
		$entry2 .= $_;
		for(my $i = 1; $i <= $cov; $i++){
			$entry1 =~ m/(transcript\/\d+)\s/;
			print STDOUT "@".$1."_$i\n";
			print STDOUT $entry2;
		}
		$line = 0;
	}
}
