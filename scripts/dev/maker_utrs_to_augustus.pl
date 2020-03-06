#!/usr/bin/perl

# script parses maker gff3 and extracts UTRS, reformats to gtf for training AUGUSTUS
# Katharina J. Hoff
# June 22nd 2018

my $id;
while(<STDIN>){
	if(m/\tmRNA\t.*Name=(g\d+\.t\d+);/){
		$id = $1;
	}elsif(m/\tfive_prime_UTR\t/){
		my @t = split(/\t/);
		print "$t[0]\t$t[1]\t5'-UTR\t$t[3]\t$t[4]\t$t[5]\t$t[6]\t$t[7]\ttranscript_id \"$id\"; gene_id \"$id\";\n";
	}elsif(m/\tthree_prime_UTR\t/){
		my @t = split(/\t/);
		print "$t[0]\t$t[1]\t3'-UTR\t$t[3]\t$t[4]\t$t[5]\t$t[6]\t$t[7]\ttranscript_id \"$id\"; gene_id \"$id\";\n";
	}
}