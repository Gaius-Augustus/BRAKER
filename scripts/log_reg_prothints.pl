#!/usr/bin/env perl

# Katharina J. Hoff
# Parse ProtHint prothint.gff (file for GeneMark-EX) to produce two classes of hints for AUGUSTUS
# Edit by Tomas Bruna: keep the al_score information in the output

use Getopt::Long;
use List::Util sum;
use POSIX;
use Math::Utils;

my $usage = <<'ENDUSAGE';

log_reg_prothints.pl     converting prothint.gff to prothint_augustus.gff

SYNOPSIS

log_reg_prothints.pl [OPTIONS] --prothint=prothint.gff --out=out.gff

OPTIONS

  --coef0 FLOAT    Default -4.00529

  --coef1 FLOAT    Default 4.73909 for mult_norm

  --coef2 FLOAT    Default 9.09026 for al_score

ENDUSAGE

my $in_file;
my $out_file;
my $coef0 = -4.00529;
my $coef1 = 4.73909;
my $coef2 = 9.09026;
my $help;

GetOptions(
    'prothint=s' => \$in_file,
    'out=s'      => \$out_file,
    'coef0=s'    => \$coef0,
    'coef1=s'    => \$coef1,
    'coef2=s'    => \$coef2,
    'help!'      => \$help
);

if ($help) {
    print $usage;
    exit(0);
}

my @lines;
open(IN, "<", $in_file) or die("ERROR: in file " . __FILE__ ." at line "
    . __LINE__ ."\nFailed to open file $in_file!\n");
@lines = <IN>;
close(IN) or die("ERROR: in file " . __FILE__ ." at line "
    . __LINE__ ."\nFailed to close file $in_file!\n");

my @mults;
foreach(@lines){
	my @t = split(/\t/);
	push(@mults, $t[5]);
}

# compute median
sub median {
  sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2;
}

my $mult_median = median(@mults);

my @normalized_mults;
foreach(@mults){
	push(@normalized_mults, ($_/$mult_median));
}

print(scalar(@normalized_mults));

open(OUT, ">", $out_file) or die("ERROR: in file " . __FILE__ ." at line "
    . __LINE__ ."\nFailed to open file $out_file for writing!\n");

for (my $i = 0; $i < scalar(@normalized_mults); $i++){
	my @t = split(/\t/, $lines[$i]);
	$t[8] =~ m/al_score=([^;]+);/;
	my $al_score = $1;
	my $y = -4.00529 + 4.73909 * $normalized_mults[$i] + 9.09026 * $al_score;
	my $class_label = 0;
	if($y >= 0.85) {
		$class_label = 2;
	}
	$t[2] =~ s/Intron/intron/;
	$t[2] =~ s/start_codon/start/;
	$t[2] =~ s/stop_codon/stop/;
	print OUT $t[0]."\t".$t[1]."\t".$t[2]."\t".$t[3]."\t".$t[4]."\t".$class_label."\t".$t[6]."\t".$t[7]."\tsrc=P;mult=$t[5];pri=4;al_score=$al_score;\n";
}

close(OUT) or die("ERROR: in file " . __FILE__ ." at line "
    . __LINE__ ."\nFailed to close file $out_file!\n");