#!/usr/bin/perl

####################################################################################################
#                                                                                                  #
# blastEx.pl - blast AUGUSTUS protein sequences based on different extrinsic files                 #
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
# | first outline                   | Simone Lange   |31.03.2015 |
# | added more sub routines         | Simone Lange   |01.04.2015 |
# | minor adjustments               | Simone Lange   |07.04.2015 |
# ----------------------------------------------------------------
 
use Getopt::Long;
use Cwd;
use File::Spec::Functions qw(rel2abs);
use Bio::Tools::Run::RemoteBlast;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SeqIO;
use File::Basename qw(dirname basename);
BEGIN{
  $0 = rel2abs($0);
  our $directory = dirname($0);
} 
use lib $directory;
use helpMod qw(find);
use Data::Dumper;

use strict;
use warnings;

my $usage = <<'ENDUSAGE';

blastEx.pl     ...

SYNOPSIS

blastEx.pl [OPTIONS] --files=augustus1.gff,augustus2.gff

  --files                     list of AUGUSTUS prediction files in gff format with protein sequence  
  
    
    
OPTIONS

    --help                          Print this help message
    --DB=path/to/blast/DB           path to local blast database
    --dir=path/to/dir               working directory
    --log=lastEx.log                log file
    --out=output.out                contains file index of highest scoring file
    --introns=hintsfile.gff         gff files with intron hints

    

Format:
  seqname <TAB> source <TAB> feature <TAB> start <TAB> end <TAB> score <TAB> strand <TAB> frame <TAB> gene_id value <TAB> transcript_id value
                           

DESCRIPTION
      
  Example:

    blastEx.pl [OPTIONS] --files=augustus1.gff,augustus2.gff

ENDUSAGE



my ($introns, $dir, @pred_files, $log, $DB, $out, $help); # necessary input files and other options
my $alpha = 1;              # factor for score (RNA-seq supported introns)
my $average_gene_length =0; # for length determination
my $average_nr_introns = 0; # average number of introns (only complete genes)
my $beta = 1;               # factor for score (bases with protein support)
my $blastfactory;
my $bool_complete = "false";# true if gene is complete, i.e. genes with start and stop codon
my $bool_good = "true";     # if true gene is good, if false, gene is bad
my $bool_intron = "false";  # true, if currently between exons, false otherwise
my $count_CDS = 0;          # count the number of CDS per gene 
my $gene_start;             # for length determination
my $good_mults = 0;         # all supported intron 'mult' entries summed up (all introns supported in gene)
my @gtf_files;              # gtf files from gff input files   
my $hintsfile;              # hints from bam2hints
my $highest_index;          # file index of the highest score
my $highest_score = '-inf'; # highest score 
my @ID;                     # gene ID
my $intron_mults = 0;       # all intron 'mult' entries summed up
my %introns;                # Hash of arrays of hashes. Contains information from intron file input. 
                            # Format $intron{seqname}{strand}[index]->{start} and ->{end}
my $length = 0;             # for length determination
my $limit = 1000;
my $mults = 0;              # current 'mult' entries
my $nr_of_bad = 0;          # number of bad genes
my $nr_of_complete = 0;     # number of complete genes, i.e. genes with start and stop codon
my $nr_of_genes = 0;        # number of all genes
my $nr_of_good = 0;         # number of good genes
my $nr_of_introns = 0;      # number of all introns in hints file
my $one_exon_gene_count = 0;# counts the number of genes which only consist of one exon
my @parameters = (-program    => "blastp",
                  -database   => "nr",
                  -expect     => '1e-5',
                  -readmethod => "SearchIO"); # parameters for blast
my @prot_files;             # protein sequence files '.aa'
my %prot_seqs;              # hash of protein sequences
my $RNAseq_sup_introns = 0; # introns supported by RNA-seq data (sum of 'mult' entries)
my %scores;                 # score parts for each input file
my $start_ID = "";          # ID of current start codon
my $sup_mults = 0;          # sum of all supported 'mult' entries (all introns, but not necessarily all introns in gene are supported)
my $true_count = 0;         # counts the number of supported CDS per gene


if(@ARGV==0){
  print "$usage\n"; 
  exit(0);
}

GetOptions( 'introns=s' => \$introns,
            'files=s'   => \@pred_files,
            'dir=s'     => \$dir,
            'out=s'     => \$out,
            'DB=s'      => \$DB,
            'log=s'     => \$log,  
            'help!'     => \$help);

if($help){
  print $usage;
  exit(0);
}

if(!defined($dir)){
 $dir = cwd(); 
}

if(!defined($out)){
 $out = "$dir/highest.index.out"; 
}

my $last_char = substr($dir, -1);
if($last_char eq "\/"){
   chop($dir);
}

if(!defined($log)){
  $log = "$dir/blastEx.log"; 
}
open (LOG, ">>".$log) or die "Cannot open file: $log\n";
  
if(defined($DB)){
  if(! -d $DB){
    print STDOUT "WARNING: Path to local BLAST database $DB does not exist. Programme will try to connect to server.\n";
    $blastfactory = Bio::Tools::Run::RemoteBlast->new(@parameters);
  }else{
    $blastfactory = Bio::Tools::Run::StandAloneBlast->new(@parameters);
  }
}else{  
  $blastfactory = Bio::Tools::Run::RemoteBlast->new(@parameters);
}

# check whether hints file is specified
if(defined($introns)){
  # check whether hints file exists
  if(! -e $introns){
    print STDOUT "WARNING: Hints file $introns does not exist. Please check.\nProgramme will only blast protein sequence files.\n";
  }else{
    $introns = rel2abs($introns);
    introns();
  }
}

# check whether prediction files exist
if(@pred_files){
  @pred_files = split(/[\s,]/, join(',',@pred_files));
  for(my $i=0; $i<scalar(@pred_files); $i++){
    if(! -e $pred_files[$i]){
      print STDERR "ERROR: Prediction file $pred_files[$i] does not exist. Please check.\n";
      exit(1);
    }
    $pred_files[$i] = rel2abs($pred_files[$i]);
    # make '*.aa' files
    getAnnoFasta($pred_files[$i]);
    if(defined($introns)){
      # make '*.gtf' files
      make_gtf($pred_files[$i]);
    }
  }
}else{
  print STDERR "ERROR: No prediction file assigned. Please check.\n";
  exit(1);
}

# get protein sequences from files and store in hash
get_prots();
# compare introns from '*.gtf' files to introns in hints file 
for(my $i=0; $i<scalar(@gtf_files); $i++){
  compare_introns($gtf_files[$i]);
  $scores{$gtf_files[$i]}{"sup_introns"} = $sup_mults;
  $scores{$gtf_files[$i]}{"pos"} = 0;
  $scores{$gtf_files[$i]}{"length"} = 0;
}

blast();

for(my $i=0; $i<scalar(@gtf_files); $i++){
  my $score = $alpha * $scores{$gtf_files[$i]}{"sup_introns"} - $intron_mults + $beta * $scores{$gtf_files[$i]}{"pos"} - $scores{$gtf_files[$i]}{"length"};
  if($score >= $highest_score){
    $highest_score = $score;
    $highest_index = $i;
  }
}

open (HIGHEST, ">$out") or die "Cannot open file: $out\n";
print HIGHEST "$highest_index";
close (HIGHEST) or die("Could not close file $out!\n");
close(LOG) or die("Could not close file $log!\n");


                           ############### sub functions ##############

# get protein sequences from files
sub get_prots{
  my $seqfile; 
  for(my $i=0; $i<scalar(@prot_files); $i++){
    print LOG "\# ".localtime.": get protein sequences from $prot_files[$i]\n";
    $seqfile = Bio::SeqIO -> new (-file => "<$prot_files[$i]",
                                  -format => "fasta");
    while(my $seq = $seqfile -> next_seq){
      if(defined($prot_seqs{$seq -> seq})){
        $prot_seqs{$seq -> seq}{"indices"} .= ",$i";
      }else{
        $prot_seqs{$seq -> seq}{"indices"} .= "$i";
        $prot_seqs{$seq -> seq}{"sobj"} = $seq;
      }
    }
  }
}

# blast sequences 
sub blast{
  foreach my $seq (keys %prot_seqs){
    print LOG "\# ".localtime.": blast ".$prot_seqs{$seq}{"sobj"} -> id."\n";
    my @file_index = split(/,/, $prot_seqs{$seq}{"indices"});
    if(scalar(@file_index) != scalar(@prot_files) || scalar(@prot_files) == 1){
      my $r = $blastfactory -> submit_blast($prot_seqs{$seq}{"sobj"});
    }
  }
  while(my @rids = $blastfactory -> each_rid){
    foreach my $rid (@rids){
      my $rc = $blastfactory -> retrieve_blast($rid); # Bio::SearchIO::blast
      if(!ref($rc)){
        if($rc < 0){
          $blastfactory -> remove_rid($rid);
        }
        print STDERR ".";
        sleep 5;
      }else{
        my $result = $rc -> next_result(); # Bio::Search::Result::BlastResult
        print LOG "\# ".localtime.": get results for ".$result -> query_name()."\n";
        $blastfactory -> remove_rid($rid);
        my $hit = $result -> next_hit;
        my $hsp = $hit -> next_hsp;
        if(defined($prot_seqs{$hsp -> query_string})){
          my @file_index = split(/,/, $prot_seqs{$hsp -> query_string}{"indices"});
          for(my $i = 0; $i<scalar(@file_index); $i++){
            $scores{$file_index[$i]}{"pos"} += $hsp -> num_conserved;
            $scores{$file_index[$i]}{"length"} += $prot_seqs{$hsp -> query_string}{"sobj"} -> length;         
          }
        }
      }
    }
  }
}

# read in introns
sub introns{
  open (INTRONS, $introns) or die "Cannot open file: $introns\n";
  print LOG "\# ".localtime.": read in introns from $introns\n";
  while(<INTRONS>){
    chomp;
    my @line = split(/\t/, $_);
    if(scalar(@line) == 9){
      $introns{$line[0]}{$line[6]}{$line[3]}{$line[4]} = $line[5];
      $nr_of_introns++;
      $intron_mults += $line[5];
    }
  }  
  close(INTRONS) or die("Could not close file $introns!\n");
}

# call getAnnoFasta.pl
sub getAnnoFasta{
  my $AUG_pred = shift;
  @_ = split(/\//, $AUG_pred);
  my $name_base = substr($_[-1],0,-4);
  my $string = find("getAnnoFasta.pl");
  my $errorfile = "$dir/errors/getAnnoFasta.$name_base.stderr";
  my $stdoutfile = "$dir/getAnnoFasta.$name_base.stdout";
  my $perlCmdString = "perl $string $AUG_pred 1>$stdoutfile 2>$errorfile";
  print LOG "\# ".localtime.": make a fasta file with protein sequences for $AUG_pred\n";
  print LOG "$perlCmdString\n\n";
  system("$perlCmdString")==0 or die("failed to execute: $!\n");
  my $prot_file = substr($AUG_pred,0,-4).".aa";
  push(@prot_files, $prot_file); 
}

# make gtf file
sub make_gtf{
  my $AUG_pred = shift;
  @_ = split(/\//, $AUG_pred);
  my $name_base = substr($_[-1],0,-4);
  my $gtf_file = substr($AUG_pred,0,-4).".gtf";
  my $errorfile = "$dir/errors/gtf2gff.$name_base.stderr";
  my $perlstring = find("gtf2gff.pl");
  my $cmdString = "cat $AUG_pred | perl -ne 'if(m/\\tAUGUSTUS\\t/){print \$_;}' | perl $perlstring --printExon --out=$gtf_file 2>$errorfile";
  print "$cmdString\n\n";
  print LOG "\# ".localtime.": make a gtf file from $AUG_pred\n";
  print LOG "$cmdString\n\n";
  system("$cmdString")==0 or die("failed to execute: $!\n");
  push(@gtf_files, $gtf_file); 
}

# compare input files to introns from hints file
sub compare_introns{
  my $gtf_file = shift;
  my $exon;                   # current exon
  my $intron_start;
  my $intron_end;
  my $prev_ID = "no_ID";

  $sup_mults = 0;
  $good_mults = 0;
  $average_nr_introns = 0;
  $nr_of_genes = 0;
  $nr_of_complete = 0;
  $length = 0;
  $count_CDS = 0;

  open (GTF, "<".$gtf_file) or die "Cannot open file: $gtf_file\n";
  while(<GTF>){
    chomp;
    my @line = split(/\t/, $_);
    if(scalar(@line) == 9){
      if(($line[2] ne "start_codon") && ($line[2] ne "stop_codon") && ($line[2] ne "CDS")){
        next;
      }else{
        @ID = split(/\s/,$line[8]); # gene_id "gene:SPAC212.11"; transcript_id "transcript:SPAC212.11.1";
        # new gene starts
        if($prev_ID ne $ID[1]){
          if($count_CDS != 0){
            $nr_of_genes++;
            reset_param();
          }
        }
        if( ($line[2] eq "start_codon" && $line[6] eq "+") || ($line[2] eq "stop_codon" && $line[6] eq "-") ){
          $gene_start = $line[3];
        # gene ends
        }elsif(($line[2] eq "stop_codon" && $line[6] eq "+") || ($line[2] eq "start_codon" && $line[6] eq "-") ){
          if($start_ID eq $ID[1]){
            $length += $line[4] - $gene_start;
            $nr_of_complete++;
            $bool_complete = "true";
          }
        # exons, CDS usw., i.e. no start or stop codon
        }elsif($line[2] eq "CDS"){
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
            $intron_start = $line[4]+1;
          }
          $count_CDS++;
        }
        $prev_ID = $ID[1];
      }
    }
  }
  $nr_of_genes++;
  reset_param(); 

  close (GTF) or die("Could not close file $gtf_file!\n");
  $average_nr_introns = $average_nr_introns / $nr_of_genes;
}

# reset parameters
sub reset_param{
  if( ($true_count + 1 ) != $count_CDS && $count_CDS != 1){
    $bool_good = "false";
  }
  if($count_CDS == 1){
    $one_exon_gene_count++;
  }

  if($bool_complete eq "true"){
    $average_nr_introns += $count_CDS - 1;
  }
  # all exons in intron file
  if($bool_good eq "true"){
    $nr_of_good++;
    $good_mults += $mults;
   # not all exons in intron file or gene incomplete
   }else{
    $nr_of_bad++;
  }
  $sup_mults += $mults;
  $count_CDS = 0;
  $true_count = 0;
  $bool_intron = "false";
  $bool_good = "true";
  $bool_complete = "false";
  $mults = 0;
}
