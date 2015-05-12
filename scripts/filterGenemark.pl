#!/usr/bin/perl

####################################################################################################
#                                                                                                  #
# filterGenemark.pl - reformats and filters the GeneMark-ET output for usage with braker.pl:       #
#                     adds double quotes around ID to match gtf format if necessary                #
#                     filters GeneMark-ET output into good and bad genes, i.e.                     #
#                     genes included and not included in introns file respectively                 #
#                                                                                                  #
# Author: Simone Lange                                                                             #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# Release date: January 7th 2015                                                                   #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
#                                                                                                  #
#################################################################################################### 


# ----------------------------------------------------------------
# | first outline from old version  | Simone Lange   |29.07.2014 |
# | changed to adapt to             |                |03.10.2014 |  
# | GeneMark-ET changes up to       |                |08.10.2014 |
# | Version 4.17 & 4.21             |                |12.01.2015 |
# | minor corrections and           |                |07.10.2014 |
# | simplifications                 |                |08.10.2014 |
# | removed '"' in 'Number of       |                |25.01.2015 |
# | genes: "nr' output              |                |           |
# | added more information to       |                |05.03.2015 |
# | standard output                 |                |           |
# | added suppress option           |                |26.03.2015 |
# | (only standard output)          |                |           |
# | more if-conditions (no division |                |           |
# | by zero errors)                 |                |           |
# ----------------------------------------------------------------
 
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(ceil);




my $usage = <<'ENDUSAGE';

filterGenemark.pl     filter GeneMark-ET files and search for "good" genes

SYNOPSIS

filterGenemark.pl [OPTIONS] genemark.gtf introns.gff

  genemark.gtf         file in gtf format
  introns.gff          corresponding introns file in gff format
  
    
    
OPTIONS

    --help                          Print this help message
    --introns=introns.gff           Corresponding intron file in gff format
    --genemark=genemark.gtf         File in gtf format
    --output=newfile.gtf            Specifies output file name. Default is 'genemark-input_file_name.c.gtf' 
                                    and 'genemark-input_file_name.f.good.gtf'
                                    and 'genemark-input_file_name.f.bad.gtf' for filtered genes included and not included in intron file respectively

Format:
  seqname <TAB> source <TAB> feature <TAB> start <TAB> end <TAB> score <TAB> strand <TAB> frame <TAB> gene_id value <TAB> transcript_id value
                           

DESCRIPTION
      
  Example:

    filterGenemark.pl [OPTIONS] --genemark=genemark.gtf --introns=introns.gff

ENDUSAGE


my ($genemark, $introns, $output_file, $help);
my $average_gene_length =0; # for length determination
my $average_nr_introns = 0; # average number of introns (only complete genes)
my $bool_complete = "false";# true if gene is complete, i.e. genes with start and stop codon
my $bool_good = "true";     # if true gene is good, if false, gene is bad
my $bool_intron = "false";  # true, if currently between exons, false otherwise
my @CDS;                    # contains coding region lines
my $file_name;              # file name
my $gene_start;             # for length determination
my $good_mults = 0;         # all supported intron 'mult' entries summed up
my $ID_new;                 # new ID with doublequotes for each gene
my @ID_old;                 # old ID without quotes
my %introns;                # Hash of arrays of hashes. Contains information from intron file input. 
                            # Format $intron{seqname}{strand}[index]->{start} and ->{end}
my $intron_mults = 0;       # all intron 'mult' entries summed up
my $length = 0;             # for length determination
my @line;                   # each input file line
my $mults = 0;              # current 'mult' entries
my $nr_of_bad = 0;          # number of bad genes
my $nr_of_complete = 0;     # number of complete genes, i.e. genes with start and stop codon
my $nr_of_genes = 0;        # number of all genes
my $nr_of_good = 0;         # number of good genes
my $one_exon_gene_count = 0;# counts the number of genes which only consist of one exon
my $start_codon = "";       # contains current start codon for "+" strand (stop codon for "-" strand)
my $start_ID = "";          # ID of current start codon (only important for files without 'stop_codon' entries)
my $stop_codon = "";        # contains current stop codon for "+" strand (start codon for "-" strand)
my $suppress;               # suppress file output
my $true_count = 0;         # counts the number of true cases per gene



if(@ARGV==0){
  print "$usage\n"; 
  exit(0);
}

GetOptions( 'genemark=s' => \$genemark,
            'introns=s'  => \$introns,
            'output=s'   => \$output_file,
            'suppress!'  => \$suppress,
            'help!'      => \$help);

if($help){
  print $usage;
  exit(0);
}

# set $genemark
if(!defined($genemark)){
  $genemark = $ARGV[0];
}

# set $introns
if(!defined($introns)){
  $introns = $ARGV[1];
}

# check whether the genemark file exists
if(!defined($genemark)){
  print STDERR "No genemark file specified. Please set a file with --genemark=genemark-ET.gtf.\n";
  exit(1),
}else{
  if(! -f "$genemark"){
    print STDERR "Genemark file $genemark does not exist. Please check.\n";
    exit(1);
  }
}

# set $file_name
if(!defined($output_file)){
  $file_name = substr($genemark,0,-4);
}else{
  $file_name = substr($output_file,0,-4);
}

# check for option
if(defined($introns)){
  # check whether the intron file exists
  if(! -f "$introns"){
    print STDERR "Intron file $introns does not exist. Please check.\n";
    exit(1);
  }else{
    introns();
  }
}
convert_and_filter();

if($nr_of_complete > 0){
  $average_gene_length = ceil($length / $nr_of_complete);
  print STDOUT "Average gene length: $average_gene_length\n";
  $average_nr_introns = $average_nr_introns / $nr_of_complete;
  print STDOUT "Average number of introns: $average_nr_introns\n";
  my $rate_genes_g = $nr_of_good / $nr_of_complete;
  print STDOUT "Good gene rate: $rate_genes_g\n";
}else{
  print STDERR "Average gene length cannot be computed since all genes are incomplete.\n";
  print STDERR "Average number of introns cannot be computed since all genes are incomplete.\n";
  print STDERR "Rate of genes cannot be computed since all genes are incomplete.\n";
}
print STDOUT "Number of genes: $nr_of_genes\n";
print STDOUT "Number of complete genes: $nr_of_complete\n";
print STDOUT "Number of good genes: $nr_of_good\n";
print STDOUT "Number of one-exon-genes: $one_exon_gene_count\n";
print STDOUT "Number of bad genes: $nr_of_bad\n";
if($intron_mults > 0){
  my $rate_mult = $good_mults / $intron_mults;
  print STDOUT "Good intron rate: $rate_mult\n";
}else{
  print STDERR "Rate of good introns cannot be computed since there are no 'mult' entries.\n";
}
if($nr_of_good > 0){ 
  my $onex_rate_g = $one_exon_gene_count / $nr_of_good;
  print STDOUT "One exon gene rate (of good genes): $onex_rate_g\n";
}else{
  print STDERR "Rate of one exon genes (of good genes) cannot be computed since only complete genes can be 'goood'.\n";
}
if($nr_of_genes > 0){
  my $onex_rate_a = $one_exon_gene_count / $nr_of_genes;
  print STDOUT "One exon gene rate (of all genes): $onex_rate_a\n";
}else{
  print STDERR "Rate of one exon genes (of all genes) cannot be computed because there are no genes in the input file with the correct format.\n";
}

if(!defined($suppress)){
  open (GENELENGTH, ">$file_name.average_gene_length.out") or die "Cannot open file: $file_name.average_gene_length.out\n";
  print GENELENGTH "$average_gene_length\t$average_nr_introns\n";
  close (GENELENGTH) or die("Could not close file $file_name.average_gene_length.out!\n");
}


                           ############### sub functions ##############

# read in introns
sub introns{
  open (INTRONS, $introns) or die "Cannot open file: $introns\n";
  while(<INTRONS>){
    chomp;
    @line = split(/\t/, $_);
    if(scalar(@line) == 9){
      $introns{$line[0]}{$line[6]}{$line[3]}{$line[4]} = $line[5];
      $intron_mults += $line[5];
    }
  }  
  close (INTRONS) or die("Could not close file $introns!\n");
}



# convert genemark file into regular gtf file with double quotes around IDs
# and split genes into good and bad ones
sub convert_and_filter{
  my $exon;                   # current exon
  my $intron_start;
  my $intron_end;
  my $prev_ID = "no_ID";
  my $output_file_good;
  my $output_file_bad;
  if(!defined($output_file)){
    $output_file = "$file_name.c.gtf";
  }
  
  if(defined($introns) && !defined($suppress)){
    $output_file_good = "$file_name.f.good.gtf";
    $output_file_bad  = "$file_name.f.bad.gtf";

    open (GOOD, ">".$output_file_good) or die "Cannot open file: $output_file_good\n";
    open (BAD, ">".$output_file_bad) or die "Cannot open file: $output_file_bad\n";
  }
  if(!defined($suppress)){
    open (OUTPUT, ">".$output_file) or die "Cannot open file: $output_file\n";
  }
  open (GENEMARK, "<".$genemark) or die "Cannot open file: $genemark\n";

  while(<GENEMARK>){
    chomp;
    @line = split(/\t/, $_);
    if(scalar(@line) == 9){
      @ID_old = split(/\s/,$line[8]);
      chop($ID_old[1]);
      chop($ID_old[3]);
      if($ID_old[1] =~m/^"\w+"$/ && $ID_old[3] =~m/^"\w+"$/){
        $ID_new = $line[8];
      }else{
        $ID_new = "$ID_old[0] \"$ID_old[1]\"\; $ID_old[2] \"$ID_old[3]\"\;";
      }
      # new gene starts
      if($prev_ID ne $ID_old[1] && defined($introns)){
        if(@CDS){
          $nr_of_genes++;
          print_gene();
        }
      }
      if( ($line[2] eq "start_codon" && $line[6] eq "+") || ($line[2] eq "stop_codon" && $line[6] eq "-") ){
        $start_codon = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$ID_new\n";
        $gene_start = $line[3];
        $start_ID = $ID_old[1];
      # gene ends
      }elsif(($line[2] eq "stop_codon" && $line[6] eq "+") || ($line[2] eq "start_codon" && $line[6] eq "-") ){
        if($start_ID eq $ID_old[1]){
          $length += $line[4] - $gene_start;
          $nr_of_complete++;
          $bool_complete = "true";
        }
        $stop_codon = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$ID_new\n";
      # exons, CDS usw., i.e. no start or stop codon
      }elsif($line[2] ne "start_codon" && $line[2] ne "stop_codon" && defined($introns)){
        if($line[2] eq "CDS"){
          if($bool_intron eq "false"){
            $intron_start = $line[4] + 1;
            $bool_intron = "true";
          }else{
            $intron_end = $line[3] - 1;
          
            # check if exons are defined in intron hash made of intron input
            if(defined($introns{$line[0]}{$line[6]}{$intron_start}{$intron_end})){
              $true_count++;
              $mults += $introns{$line[0]}{$line[6]}{$intron_start}{$intron_end};
            }
            $intron_start = $line[4] + 1;
          }
          $exon = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$ID_new";
          push(@CDS, $exon);
          $exon = "$line[0]\t$line[1]\texon\t$line[3]\t$line[4]\t0\t$line[6]\t.\t$ID_new";
          push(@CDS, $exon);
        }
      }
      if(!defined($suppress)){
        print OUTPUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$ID_new\n";
      }
      $prev_ID = $ID_old[1];
    }
  }
  if(defined($introns)){
    $nr_of_genes++;
    print_gene(); # print last gene, since print_gene() was only executed after the ID changed
    if(!defined($suppress)){
      close (BAD) or die("Could not close file $output_file_bad!\n");
      close (GOOD) or die("Could not close file $output_file_good!\n");
    }
  }

  close (GENEMARK) or die("Could not close file $genemark!\n");
  if(!defined($suppress)){
    close (OUTPUT) or die("Could not close file $output_file!\n");
  }
}


# print genes in corresponding file (good or bad)
sub print_gene{
  my $size = scalar(@CDS) / 2; 
  if( ($true_count + 1 ) != $size && $size != 1){
    $bool_good = "false";
  }
  if($size == 1){
    $one_exon_gene_count++;
  }
  if($bool_complete eq "true"){
    $average_nr_introns += $size - 1;
  }
  # all exons in intron file
  if($bool_good eq "true" && $bool_complete eq "true"){
    $nr_of_good++;
    $good_mults += $mults;
    if(!defined($suppress)){
      print GOOD "$start_codon";
      foreach (@CDS){
        print GOOD "$_\n";
      }
      print GOOD "$stop_codon";
    }
   # not all exons in intron file or gene incomplete
   }else{
    $nr_of_bad++;
    if(!defined($suppress)){
      print BAD "$start_codon";
      foreach (@CDS){
        print BAD "$_\n";
      }
      print BAD "$stop_codon";
    }
  }
  @CDS =();
  $true_count = 0;
  $start_codon = "";
  $stop_codon = "";
  $bool_intron = "false";
  $bool_good = "true";
  $bool_complete = "false";
  $mults = 0;
}


