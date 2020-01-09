#!/usr/bin/env perl

####################################################################################################
#                                                                                                  #
# filterIntronsFindStrand.pl - finds corresponding strand for introns in fasta file                #
#                              optionally set the score column to the 'mult' entry with --score    #
#                                                                                                  #
# Author: Simone Lange & Katharina J. Hoff                                                         #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# Release date: December 11th 2019                                                                 #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
#                                                                                                  #
####################################################################################################

# ------------------------------------------------------------------
# | file creation and findStrand()    | Simone Lange   |06.10.2014 |
# | add getScore() for score option   |                |07.10.2014 |
# | add error message if sequence     |                |23.01.2015 |
# | name of hints and fasta file do   |                |           |
# | not match -> program   stops then |                |           |
# | excempt already stranded hints    | Katharina Hoff |11.12.2019 |
# | output both strands if ambigous   |                |           |
# ------------------------------------------------------------------

use strict;
use warnings;
use Getopt::Long;



my $usage = <<'ENDUSAGE';

filterIntronsFindStrand.pl     find corresponding strand for introns from two input files genome.fa and introns.gff

SYNOPSIS

filterIntronsFindStrand.pl genome.fa introns.gff [OPTIONS] > introns.s.f.gff

  genome.fa           DNA file in fasta format
  introns.gff         corresponding introns file in gff format


OPTIONS

    --help                          Print this help message
    --allowed=gtag,gcag,atac        Allowed acceptor and donor splice site types
    --score                         Set score to 'mult' entry or '1', if the last column does not contain a 'mult' entry
    --genome=genome.fa              see above
    --introns=introns.gff           see above




DESCRIPTION

  Example:

    filterIntronsFindStrand.pl genome.fa introns.gff [OPTIONS] > introns.s.f.gff

ENDUSAGE


my ($genome, $introns, @allowed, $mult_score, $help);
my %annos; # keys: sequences, elements: annotations
my $seqname;
my $seq;

if(@ARGV==0){
  print "$usage\n";
  exit(0);
}

GetOptions( 'introns=s' => \$introns,
            'genome=s'  => \$genome,
            'score!'    => \$mult_score,
            'allowed=s' => \@allowed,
            'help!'     => \$help);

if($help){
  print $usage;
  exit(0);
}

# set $genome
if(!defined($genome)){
  $genome = $ARGV[0];
}

# set $introns
if(!defined($introns)){
  $introns = $ARGV[1];
}

# set allowed splice site types
if(@allowed){
  @allowed = split(/[\s,]/, join(',',@allowed));
}else{
  @allowed = ("gtag", "gcag", "atac");
}

# check whether the files exist
if(! -f "$genome"){
  print "Genome file $genome does not exist. Please check.\n";
  exit(1);
}

if(! -f "$introns"){
  print "Introns file $introns does not exist. Please check.\n";
  exit(1);
}

# genome file in fasta format
open (FASTA, "<".$genome) or die "Cannot open file: $genome\n";
while(<FASTA>) {
  chomp;
  if(m/^>(.*)/){
    if(defined($seqname)){
      $annos{$seqname} = $seq;
    }
    $seqname = $1;
    $seq = "";
  }else{
    $seq .= $_;
  }
}
$annos{$seqname} = $seq;
close(FASTA) or die("Could not close fasta file $genome!\n");

# introns hintsfile in gff format
open (INTRONS, "<".$introns) or die "Cannot open file: $introns\n";
$/="\n";
while(<INTRONS>){
    chomp;
    my @line = split(/\t/, $_);
    if(scalar(@line) == 9){
        if($line[6] =~ m/\./){
            my $strand = findStrand($line[0], $line[3], $line[4]);
            if( $strand eq "+" || $strand eq "-" || $strand eq "b" ) {
                my $score;
                if($mult_score){
                    $score = getScore($line[8]);
                }else{
                    $score = $line[5];
                }
                if($strand eq "+" || $strand eq "-"){
                    print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$score\t$strand\t$line[7]\t$line[8]\n";
                }elsif($strand eq "b"){
                    print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$score\t+\t$line[7]\t$line[8]\n";
                    print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$score\t-\t$line[7]\t$line[8]\n";
                }
            }
        }else{
            print $_."\n";
        }
    }
}
close(INTRONS) or die("Could not close introns file $introns!\n");


                           ############### sub functions ##############


# find strand for introns
# look up start and end position and check if it matches allowed splice site patterns
sub findStrand{
    my $seqname = shift;
    my $start = shift;
    my $end = shift;
    my $type;
    my $reverse;
    my $return_type;
    my $has_plus = 0;
    my $has_minus = 0;
    if(defined($annos{$seqname})){
        $type = lc(substr($annos{$seqname}, $start-1,2)).lc(substr($annos{$seqname}, $end-2,2));
        $reverse = reverse($type);
        $reverse =~ tr/agct/tcga/;
        foreach (@allowed){
            if($_ eq $type){
                $has_plus = 1;
            }
            if($_ eq $reverse){
                $has_minus = 1;
            }
        }
        if( ($has_minus == 1) && ($has_plus == 1) ){
            $return_type = "b";
        }elsif($has_minus == 1){
            $return_type = "-";
        }elsif($has_plus == 1){
            $return_type = "+";
        }else{
            $return_type = ".";
        }
        return $return_type;
    }else{
        print STDERR "WARNING: '$seqname' does not match any sequence in the fasta file. Maybe the two files do not belong together.\n";
        return ".";
    }
}

# get score from mult entry
sub getScore{
  my $column = shift;
  my $score;
  if($column =~ m/mult=(\d+)/){
    $score = $1;
  }else{
    $score = 1;
  }
  return $score;
}



