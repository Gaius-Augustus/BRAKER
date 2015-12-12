#!/usr/bin/perl
####################################################################################################
#                                                                                                  #
# IntronID.pl - changes last column of intron file to ID=number_mult=score                         #
#               to see the weight of each intron in the Browser                                    #
#                                                                                                  #
# Author: Simone Lange                                                                             #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# Release date: July 10th 2015                                                                     #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
#                                                                                                  #
#################################################################################################### 


# ----------------------------------------------------------------
# | creation of file                | Simone Lange   |16.12.2014 |
# ----------------------------------------------------------------
#

use Getopt::Long;
use strict;
use warnings;

my $usage = <<'ENDUSAGE';

IntronID.pl     add ID=mult to intron hints for usage with JBrowse

SYNOPSIS

IntronID.pl [OPTIONS] genemark.gtf introns.gff

  introns.gff          hints file in gff format
  

OPTIONS

    --help                          Output this help message          
    --in=introns.gff                Intron hints in gff format
    --out=new.introns.gff           Same as input plus 'D=nr_mult=score' entry in last column


Format:
  seqname <TAB> source <TAB> feature <TAB> start <TAB> end <TAB> score <TAB> strand <TAB> frame <TAB> ID=nr_mult=score
                           

DESCRIPTION
      
  Example:

    addGeneLine.pl [OPTIONS] --in=introns.gff  --out=new.introns.gff

ENDUSAGE


my $help;       # print usage 
my $i = 1;      # intron ID index
my $in;         # intron file name
my $out;        # output file name    
my $score;      # contains multiplicity as new score 

if(@ARGV==0){
  print "$usage\n"; 
  exit(0);
}

GetOptions( 'in=s'   => \$in,
      	    'out=s'  => \$out,
            'help!'  => \$help);

if($help){
  print $usage;
  exit(0);
}

# set $in
if(!defined($in)){
  $in = $ARGV[0];
}

# check whether the file exists
if(! -f "$in"){
  print "File $in does not exist. Please check.\n";
  exit(1);
}else{
  if(!defined($out)){
    print "No output file assigned. Please check.\n";
    exit(1);
  }
  open(IN, "<", $in) or die("Could not open file $in!\n");
  open(OUT, ">", $out) or die("Could not open file $out!\n");
  while(<IN>){
    chomp;
    my @line =split(/\t/,$_);
    $score = getScore($line[8]);
    print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\tID=$i\_mult=$score;\n";
    $i++;
  }

    close(IN) or die("Could not close file $in!\n");
    close(OUT) or die("Could not close file $out!\n");
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
