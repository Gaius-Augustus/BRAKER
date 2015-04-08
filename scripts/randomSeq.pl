#!/usr/bin/perl

####################################################################################################
#                                                                                                  #
# randomSeq.pl - get random sequence parts of input fasta file                                     #
#                                                                                                  #
# Author: Simone Lange                                                                             #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# Release date: April 07th 2015                                                                    #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
#                                                                                                  #
#################################################################################################### 


# ----------------------------------------------------------------
# | first outline                   | Simone Lange   |06.04.2015 |
# | some adjustments to the         | Simone Lange   |07.04.2015 |
# | subroutines                     |                |           |
# ----------------------------------------------------------------
 
use Getopt::Long;
use Cwd;
use POSIX qw(floor);
use File::Spec::Functions qw(rel2abs);
use File::Basename qw(dirname basename);

use strict;
use warnings;

my $usage = <<'ENDUSAGE';

blastEx.pl     ...

SYNOPSIS

randomSeq.pl [OPTIONS] --genome=genome.fa


  --genome=genome.fa          fasta file with DNA sequences  
  
    
    
OPTIONS

    --help                          Print this help message
    --GMET=genemark.gtf             GeneMark-ET prediction
    --dir=pathh/to/dir              working directory
    --log=lastEx.log                log file
    --out=newGenome.fa              new genome file in fasta format


                          

DESCRIPTION
      
  Example:

randomSeq.pl [OPTIONS] --genome=genome.fa 

ENDUSAGE



my ($dir, $GMET, $genome, $log, $out, $help); # necessary input files and other options
my $gene_start;             # for storing GeneMark-ET genes
my @GMET_bounds;            # start and end position of genes
my $GMET_genes = 0;         # number of genes in GeneMark-ET prediction
my @ID;                     # gene ID
my $max_length = 0;         # maximal sequence length
my $min_size = 100000;      # minimal sequence part length
my $min_genes = 200;        # minimal number of genes in sequence parts
my @new_seq;                # new fasta sequences
my $nr_of_genes = 0;        # number of genes in sequence part
my %seqs;                   # sequence information
my @seqInfo;                # sequence length and name
my $start_ID = "";          # ID of current start codon
my $whole_length = 0;       # whole length of all sequences

if(@ARGV==0){
  print "$usage\n"; 
  exit(0);
}

GetOptions( 'genome=s'  => \$genome,
            'GMET=s'    => \$GMET,
            'dir=s'     => \$dir,
            'log=s'     => \$log,
            'out=s'     => \$out, 
            'help!'     => \$help);

if($help){
  print $usage;
  exit(0);
}

if(!defined($dir)){
 $dir = cwd(); 
}

my $last_char = substr($dir, -1);
if($last_char eq "\/"){
   chop($dir);
}

if(!defined($log)){
  $log = "$dir/blastEx.log"; 
}
open (LOG, ">>".$log) or die "Cannot open file: $log\n";
  
if(!defined($out)){
  $out = "new.genome.fa"; 
}
# check whether GeneMark-ET file is specified
if(defined($GMET)){
  # check whether GeneMark-ET file exists
  if(! -e $GMET){
    print STDOUT "WARNING: GeneMark-ET file $GMET does not exist. Please check.\nProgramme will split fasta file without checking the number of genes.\n";
  }else{
    $GMET = rel2abs($GMET);
    print LOG "\# ".localtime.": Get GeneMark-ET genes from $GMET\n";
    GMET();
  }
}

# check whether genome file exists
if(defined($genome)){
  # check whether genome file exists
  if(! -e $genome){
    print STDERR "ERROR: Genome file $genome does not exist. Please check.\n";
    exit(1);
  }
  $genome = rel2abs($genome);

  if($GMET_genes > $min_genes){
    print LOG "\# ".localtime.": Get fasta sequences from $genome\n";
    get_fasta();
    print LOG "\# ".localtime.": Get fasta sequence parts\n";
    while($nr_of_genes <= $min_genes){
      get_random_seq();
    }
    print LOG "\# ".localtime.": Print fasta sequence parts to $out\n";
    open (OUT, ">".$out) or die "Cannot open file: $out\n";
    for(my $i=0; $i<scalar(@new_seq); $i++){
      my $name = $i + 1;
      print OUT ">$name\n";
      while(length($new_seq[$i]) > 50 ){
        print OUT substr($new_seq[$i],0 ,50)."\n";
        substr($new_seq[$i],0, 50, "");
      }
      print OUT $new_seq[$i]."\n";
    }
    close(OUT) or die("Could not close file $out!\n");
  }else{
    print STDOUT "WARNING: Number of genes in $GMET is less ($GMET_genes) than $min_genes. Programme will use the whole fasta file $genome\n";
    print LOG "\# ".localtime.": copy fasta file\n";
    print LOG "cp $genome $out\n\n";
    my $cmdString = "cp $genome $out";
    system("$cmdString")==0 or die("failed to execute: $!\n");
  }
  close(LOG) or die("Could not close file $log!\n");
}else{
  print STDERR "ERROR: No genome file assigned. Please check.\n";
  exit(1);
}



                           ############### sub functions ##############

# read in GeneMark-ET genes
sub GMET{
  my $prev_ID = "no_ID";
  open (GMET, $GMET) or die "Cannot open file: $GMET\n";
  print LOG "\# ".localtime.": read in genes from $GMET\n";
  while(<GMET>){
    chomp;
    my @line = split(/\t/, $_);
    if(scalar(@line) == 9){
      if(($line[2] ne "start_codon") && ($line[2] ne "stop_codon")){
        next;
      }else{
        @ID = split(/\s/,$line[8]);
        # new gene starts
        if( ($line[2] eq "start_codon" && $line[6] eq "+") || ($line[2] eq "stop_codon" && $line[6] eq "-") ){
          $gene_start = $line[3];
          $start_ID = $ID[1];
        # gene ends
        }elsif(($line[2] eq "stop_codon" && $line[6] eq "+") || ($line[2] eq "start_codon" && $line[6] eq "-") ){
          if($start_ID eq $ID[1]){
            $GMET_genes++;
            push(@{$seqs{$line[0]}}, $gene_start, $line[4]);
          }
        }
      
        $prev_ID = $ID[1];
      }
    }
  }
  foreach my $seq (keys %seqs){
    @{$seqs{$seq}} = sort{$a  <=> $b} @{$seqs{$seq}};
  }
  close(GMET) or die("Could not close file $GMET!\n");
}


sub get_fasta{
  open (FASTA, $genome) or die "Cannot open file: $%hashgenome\n";
  print LOG "\# ".localtime.": read in DNA sequence from $genome\n";
  $/ = ">";
  while(<FASTA>){
    s/>$//;                           # see getAnnoFasta.pl
    next unless m/\S+/;               # see getAnnoFasta.pl
    /(.*)\n/;                         # see getAnnoFasta.pl
    my $seqname = $1;                 # see getAnnoFasta.pl
    my $sequencepart = $';            # see getAnnoFasta.pl
    $seqname =~ s/\s.*//;             # seqname only up to first white space (see getAnnoFasta.pl)
    $sequencepart =~ s/\s//g;         # see getAnnoFasta.pl
    my %hash = ("name"   => $seqname,
                "seq"    => $sequencepart,
                "start"  => 0);
    push(@seqInfo, \%hash);
    $whole_length += length($sequencepart);
  }
  $/ = "\n";
  close(FASTA) or die("Could not close file $genome!\n");

  @seqInfo = sort{length($a->{"seq"}) <=> length($b ->{"seq"})} @seqInfo;
}


sub get_random_seq{
  my $rand = rand();
  my $low_boundary = 0;
  my $up_boundary = length($seqInfo[0] -> {"seq"}) / $whole_length;
  my $index = 0;
  for(my $i=0; $i<scalar(@seqInfo)-1; $i++){
    if($low_boundary <= $rand && $rand < $up_boundary){
      $index = $i;
      last;
    }else{      
      $low_boundary = $up_boundary;
      $up_boundary = $up_boundary + length($seqInfo[$i+1] -> {"seq"}) / $whole_length;
    }
  }
  
  my $subseq;
  my $part_seq;
  my $part_start;
  my $start = 0;
  if(length($seqInfo[$index] -> {"seq"}) >= $min_size){
    $rand = floor(rand(length($seqInfo[$index] -> {"seq"}) + 1));
    $start = $rand;
    if($rand > length($seqInfo[$index] -> {"seq"}) - $min_size){
      $subseq = substr($seqInfo[$index] -> {"seq"}, - $min_size);
      $seqInfo[$index] -> {"seq"} = substr($seqInfo[$index] -> {"seq"}, $min_size);
    }else{
      $subseq = substr($seqInfo[$index] -> {"seq"}, $rand, $min_size);
      $part_seq = substr($seqInfo[$index] -> {"seq"}, $rand + $min_size + 1);
      $part_start = $rand + $min_size + 1;
      if($rand != 0){
        $seqInfo[$index] -> {"seq"} = substr($seqInfo[$index] -> {"seq"}, 0, $rand - 1);
      }else{
        $seqInfo[$index] -> {"seq"} = "";
      }
      my %hash = ("name"   => $seqInfo[$index] -> {"name"},
                  "seq"    => $part_seq,
                  "start"  => $part_start);
      push(@seqInfo, \%hash);     
    }
  }else{
#    $subseq = $seqInfo[$index] -> {"seq"};
#    $seqInfo[$index] -> {"seq"} = "";
  }
  $whole_length -= length($subseq);

  if(length($subseq) > 0){
    if($seqInfo[$index] -> {"start"} + $start + $min_size < $seqs{$seqInfo[$index] -> {"name"}}[1]){
      $nr_of_genes += 0;
    }else{
      $low_boundary = 0;
      my $i = 0;
      
      until($seqInfo[$index] -> {"start"}  + $start < $low_boundary){
        $low_boundary = $seqs{$seqInfo[$index] -> {"name"}}[$i];
        $i++;
      }
      $i -= 1;
      if($seqInfo[$index] -> {"start"} + $start + length($subseq) < $low_boundary){
        $nr_of_genes += 0;
      }else{
        my $j = $i;

        $up_boundary = $seqs{$seqInfo[$index] -> {"name"}}[$j];
        while($seqInfo[$index] -> {"start"} + $start + length($subseq) > $up_boundary && $j < scalar(@{$seqs{$seqInfo[$index] -> {"name"}}})){
          $up_boundary = $seqs{$seqInfo[$index] -> {"name"}}[$j];
          $j++;
        }
        if($j == scalar(@{$seqs{$seqInfo[$index] -> {"name"}}})){
          $j -= 1;
        }else{
          $j -= 2;
        }
        
        if($j % 2){
          $j += 1;
        }
        if($i % 2){
          $i += 1;
        }
        $nr_of_genes += ($j - $i) / 2;
      }
    }

    push(@new_seq, $subseq);
  }
  @seqInfo = sort{length($a -> {"seq"}) <=> length($b -> {"seq"})} @seqInfo;
}

