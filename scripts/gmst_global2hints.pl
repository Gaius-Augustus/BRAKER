#!/usr/bin/env perl

# Convert GeneMark-ST global genomic genes into hints for AUGUSTUS
# Author: Katharina J. Hoff
# Date: December 1st 2021
# Contact: katharina.hoff@uni-greifswald.de
# Usage:
# gmst_global2hints.pl in.gtf > hints.out

if(not(scalar(@ARGV)==1)){
	print "Wrong usage!\n";
	print "gmst_global2hints.pl in.gtf > hints.out\n";
	exit(1);
}

my $gmstF = $ARGV[0];
my $line;
my $prev;
my $txid;

open(GMST, "<", $gmstF) or die("Could not open file $gmstF!\n");
while(<GMST>){
	$line = $_;
	@t = split(/\t/); 
	if($t[2] eq "CDS"){
 		if(not(defined($txid))){
 			$txid="";
 		} 
 		$line=~m/transcript_id "(\S+)"/;
		if(not($txid eq $1)){
			$txid=$1;
			if($t[3]+7 < $t[4]){
				$t[3] = $t[3]+7; 
				$t[2] = "CDSpart"; 
			}
			if(defined($prev)){
				@p = split(/\t/, $prev); 
				if($p[4]-7 > $p[3]){
					$p[4] = $p[4]-7;
					$p[2] = "CDSpart";
				}
				$prev = $p[0]."\t".$p[1]."\t".$p[2]."\t".$p[3]."\t".$p[4]."\t".$p[5]."\t".$p[6]."\t".$p[7]."\t".$p[8];
				if($p[3]>$p[4]){
					print("Problem: ", $prev);
				}
			}
		}
		if(defined($prev)){
			print $prev;
		}
		$prev = "$t[0]\tGMST\t$t[2]\t$t[3]\t$t[4]\t$t[5]\t$t[6]\t$t[7]\tsrc=M;pri=4;grp=$txid;\n";
	}
}
close(GMST) or die("Could not close file $gmstF!\n");
