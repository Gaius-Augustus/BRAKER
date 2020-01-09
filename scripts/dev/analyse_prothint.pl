#!/usr/bin/env perl

# Katharina Hoff
# Label prothint hints as TP/FP

# input combset file contains only TP introns list created by overlapStat.pl
# format:
# Locus\tstart\tend\tstrand\n

my $combF =  $ARGV[0]; # combset.11.lst that contains seq_start_stop_strand of true positive introns from ProtHint

my $protF = $ARGV[1]; # prothint.gff

my %true;
open(COMB, "<", $combF) or die("Could not open file $combF!\n");
while (<COMB>){
	chomp;
	$true{$_} = 1;
}
close(COMB) or die("Could not close file $combF!\n");

print("mult\tal_score\ttop_or_not\tlabel\n");

open(PROT, "<", $protF) or die("Could not open file $protF!\n");
while(<PROT>){
	if(m/Intron/){
	my @t = split(/\t/);
	my $id_string = "$t[0]_$t[3]_$t[4]_$t[6]";
	print "$t[5]\t";
	$t[8] =~ m/al_score=([^;]+);/;
	print $1."\t";
	if($t[8] =~ m/topProt=TRUE/){
		print "topProt\t";
	}else{
		print "notTop\t";
	}
	if(defined($true{$id_string})){
			print "T\n";
		}else{
			print "F\n";
		}
	}
}
close(PROT) or die("Could not close file $protF!\n");
