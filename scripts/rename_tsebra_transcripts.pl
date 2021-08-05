#!/usr/bin/env perl

# Katharina J. Hoff
# July 2nd 2021

# This script renames TSEBRA transcripts & genes with a given prefix.
#
# Call:
# cat tsebra.gtf | rename_tsebra_transcripts.pl prefix > fixed.gtf

my $stem = $ARGV[0]."_g"; # stem of gids
my %gi = ();

my $g = 0;
my $t = 0;
my $printall = 0;
while(<STDIN>){
    $_ =~ s/\r//g;
    $line = $_;
    if($_ =~ m/transcript_id \"([^"]+)\"; gene_id \"([^"]+)\"/){
	$ctxid = $1;
	$cgid = $2;
	$printall = 1;
    }elsif($_ =~ m/\ttranscript\t/){
	chomp;
	@a = split(/\t/, $_);
	$ctxid = $a[8];
	$printall = 0;
    }else{
	print("Nothing matches: " + $_ + "\n");
    }
    if(not(defined($gi{$cgid}))){
      	$gi{$cgid} = ();
	$t = 0;
      	$g++;
    }
    if(not(defined($gi{$cgid}{$ctxid}))){
	$gi{$cgid}{$ctxid} = 1;
	$t++;
    }
    if($printall == 1){
	$line =~ s/transcript_id \"[^"]+\"; gene_id \"[^"]+\";/transcript_id \"$stem$g\.t$t\"; gene_id \"$stem$g\";/;
    }else{
	@b = split(/\t/, $line);
	$line = "$b[0]\t$b[1]\t$b[2]\t$b[3]\t$b[4]\t$b[5]\t$b[6]\t$b[7]\ttranscript_id \"$stem$g\.t$t\";\n";
    }
    print $line;
}
