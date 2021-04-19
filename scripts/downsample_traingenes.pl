#!/usr/bin/env perl

# Author: Katharina J. Hoff and Mario Stanke
# Contact: katharina.hoff@uni-greifswald.de
# Last modification: August 29th 2019

# This script samples training gene gene structures according to their intron number
# Original observation is: if training gene structures originate from GeneMark-EX,
# downsampling single-exon genes and genes with few introns (1-5) improves accuracy
# of AUGUSTUS after training on those gene structures.

use strict;
use warnings;

use Getopt::Long;

my $usage = <<'ENDUSAGE';

downsample_traingenes.pl     sample training gene structures from GeneMark-EX

SYNOPSIS

downsample_traingenes.pl --in_gtf=traingenes.gtf --out_gtf=out.gtf

	traingenes.gtf   training gene structure file in gtf format (e.g. from
	                 GeneMark-EX).
	out.gtf          output file with downsampled training gene structures


OPTIONS

    --help                          Print this help message
    --version                       Print version of script
    --lambda=f                      Parameter lambda of Poisson distribution
                                    (default value is 0)
    --intron_num_lst=s              File with intron numbers per gene (selected)

DESCRIPTION

  Example:

   downsample_traingenes.pl --in_gtf=traingenes.gtf --out_gtf=out.gtf --lambda=1

ENDUSAGE

srand 10; # make downsampling reproducible

my $lambda = 0;
my $in_gtf;
my $out_gtf;
my $help;
my $version = 1.0;
my $print_version;
my $intron_num_lst;

GetOptions(
    'in_gtf=s' 		  => \$in_gtf,
    'out_gtf=s'       => \$out_gtf,
    'lambda=s'        => \$lambda,
    'help!'           => \$help,
    'version!'        => \$print_version,
    'intron_num_lst=s'=> \$intron_num_lst
    );

if(!$in_gtf) {
    die ("Input file not defined (--int_gtf=file)!\n$usage");
}

if(!$out_gtf) {
    die ("Output file not defined (--out_gtf=file)!\n$usage");
}

if ($help) {
    print $usage;
    exit;
}

if ($print_version) {
    print "Version $version\n";
    exit;
}

my %nIntrons;
my %tx;

####################### read gtf file ##########################################

open(GTF, "<", $in_gtf) or die("Could not open file $in_gtf!\n");
while(<GTF>){
    if ( $_ =~ m/transcript_id \"(\S+)\"/) {
	my $txid = $1;
	push(@{$tx{$txid}}, $_);
	if( $_ =~ m/\tCDS\t/ ) {
	    if(not(defined($nIntrons{$txid}))) {
		$nIntrons{$txid} = 0;
	    }else {
		$nIntrons{$txid}++;
	    }
	}
    }
}
close(GTF) or die("Could not close file $in_gtf!\n");

####################### compute F(X) ###########################################

# always keep genes with more than this number:
my $max_intron_number = 5;
my @F; # distribution function
for (my $i = 0; $i <= $max_intron_number; $i++ ) {
    if ( $i == 0 ) {
	$F[$i] = P_X_is_k($i, $lambda);
    }else{
	$F[$i] = $F[$i-1] + P_X_is_k($i, $lambda);
    }
}

####################### Sample genes ###########################################

if($intron_num_lst) {
	open (LST, ">", $intron_num_lst) or die ("Could not open file $intron_num_lst!\n");
}

open(OUT, ">", $out_gtf) or die("Could not open file $out_gtf!\n");

my $min_single_exon_genes = 20;
my $single_exon_gene_counter = 0;
my %intronNumPrinted;

foreach my $txid (sort(keys %nIntrons)) {
    my $intronNum = $nIntrons{$txid};
    if( $intronNum == 0 ) {
	$single_exon_gene_counter++;
    }
    my $u = rand(1);
    my $index = $intronNum;
    # For intron numbers between 0 and $max_intron_number compare a random Poisson distributed 
    # random variable N with the intron number n and keep the gene unless N<n.
    # In effect, genes with few introns are downsampled.
    if($index > $max_intron_number || $u <= $F[$index]
       || ($intronNum == 0 && $single_exon_gene_counter <= $min_single_exon_genes )){
	foreach (@{$tx{$txid}}) {
	    print OUT $_;
	    if($intron_num_lst && not(defined($intronNumPrinted{$txid}))) {
		print LST $intronNum."\t".$txid."\n";
		$intronNumPrinted{$txid} = 1;
	    }
	}
    }
}
close (OUT) or die ("Could not close file $out_gtf!\n");


if($intron_num_lst) {
    close(LST) or die ("Could not close fiel $intron_num_lst!\n");
}

if($single_exon_gene_counter < 20){
    print STDOUT "WARNING: Number of single exon training genes is smaller "
      ."than 20. It is: $single_exon_gene_counter!\n";
}
####################### P_X_is_k ###############################################
# Computes the P(X=k), currently with Poisson distribution
################################################################################

sub P_X_is_k {
    my $k = shift;
    my $L = shift;
    my $e = exp(1);
    my $result = ( $e ** (-1 * $L) ) * ( ( $L ** $k ) / ( factorial ($k) ) );
    return $result;
}

####################### factorial ##############################################
# Computes factorial(x)
################################################################################

sub factorial {
    my $n = shift;
    my $result = 1;
	for (my $i = 1; $i <= $n; $i++) {
	    $result = $result * $i;
	}
    return $result;
}
