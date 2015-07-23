#!/usr/bin/perl
####################################################################################################
#                                                                                                  #
# addGeneLine.pl - adds gene and transcript line to input files like GeneMark-ET prediction        #
#                  or reference annotation to make them compatible with JBrowse                    #
#                                                                                                  #
#                                                                                                  #
# Author: Simone Lange                                                                             #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# Release date: April 22nd 2015                                                                    #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
#                                                                                                  #
#################################################################################################### 


# ----------------------------------------------------------------
# | first outline                   | Simone Lange   |03.10.2014 |
# | changed to adapt to             |                |           |  
# | GeneMark-ET changes up to       |                |           |
# | Version 4.15                    |                |           |
# | minor corrections and           |                |07.10.2014 |
# | simplifications                 |                |           |
# ----------------------------------------------------------------
# 
use strict;
use warnings;
use Getopt::Long;
use File::Compare;
use Data::Dumper;
use POSIX qw(ceil);




my $usage = <<'ENDUSAGE';

addGeneLine.pl     adds gene and transcript line to input file

SYNOPSIS

addGeneLine.pl [OPTIONS] input.gtf

  input.gtf         file in gtf format
  

OPTIONS

    --help                          Output this help message
    --in=input.gtf                  Input file in gtf or gff3 format
    --out=newfile.gtf               Same as input plus 'gene' and 'transcript' lines      


Format:
  seqname <TAB> source <TAB> feature <TAB> start <TAB> end <TAB> score <TAB> strand <TAB> frame <TAB> gene_id value <TAB> transcript_id value
                           

DESCRIPTION
      
  Example:

    addGeneLine.pl [OPTIONS] --in=input.gtf --out=newfile.gtf

ENDUSAGE


my $currentEnd;             # current end of input line
my @entries;                # contains coding region lines
my $gene_line = "";         # additional gene line entry
my $help;                   # print usage 
my $in;                     # prediction file name
my @line;                   # each input file line
my $out;                    # Same as input plus 'gene' and 'transcript' lines
my $transcript_line = "";   # additional transcript line entry




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
  addGene();
}



                           ############### sub functions ##############


# add gene and transcript line based on first and last CDS entry of each transcript
sub addGene{
  my $prev_ID = "no_ID";

  open(IN, "<", $in) or die("Could not open file $in!\n");
  open(OUT, ">", $out) or die("Could not open file $out!\n"); 

  while(<IN>){
    chomp;
    if(length($_)>0){
      @line = split(/\t/, $_);
      next if (@line<8);
      if($line[2] eq "CDS"){
        if($prev_ID eq "no_ID"){
          $gene_line = "$line[0]\t$line[1]\tgene\t$line[3]\t";
          $transcript_line = "$line[0]\t$line[1]\ttranscript\t$line[3]\t";
        }
        if($prev_ID ne $line[8] && $prev_ID ne "no_ID"){
          $gene_line .= $currentEnd;
          $transcript_line .= $currentEnd;
          if(@entries){
            print_gene();
          }
          
          $gene_line = "$line[0]\t$line[1]\tgene\t$line[3]\t";
          $transcript_line = "$line[0]\t$line[1]\ttranscript\t$line[3]\t";
        }
        
        $prev_ID = $line[8];
        $currentEnd = "$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\n";
      }
      if($line[2] ne "5UTR" && $line[2] ne "3UTR"){ # UTR entries were not predicted under BRAKER1 -> remove them so that they do not "bother" in the genome browser 
        push(@entries, $_);
      }
    }
  }
  $gene_line .= $currentEnd;
  $transcript_line .= $currentEnd;
  print_gene(); # print last gene, since print_gene() was only executed after the ID changed

  close(IN) or die("Could not close file $in!\n");
  close(OUT) or die("Could not close file $out!\n");
}




sub print_gene{
  print OUT "$gene_line";
  print OUT "$transcript_line";
  foreach (@entries){
    print OUT "$_\n";
  }
  @entries =();
}


