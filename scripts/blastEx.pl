#!/usr/bin/perl

####################################################################################################
#                                                                                                  #
# blastEx.pl - blast AUGUSTUS protein sequences based on different extrinsic files                 #
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
# | first outline                   | Simone Lange   |31.03.2015 |
# | added more sub routines         |                |01.04.2015 |
# | minor adjustments               |                |07.04.2015 |
# | score adjustments               |                |05.05.2015 |
# | made blast forkable             |                |08.05.2015 |
# ----------------------------------------------------------------
 
use Getopt::Long;
use Cwd;
use POSIX qw(log10 ceil);
use File::Spec::Functions qw(rel2abs);
use Statistics::Basic qw(:all);
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::SeqIO;
use Parallel::ForkManager;
use File::Basename qw(dirname basename);
BEGIN{
  $0 = rel2abs($0);
  our $directory = dirname($0);
} 
use lib $directory;
use helpMod qw(find uptodate);
use Data::Dumper;

use strict;
use warnings;

my $usage = <<'ENDUSAGE';

blastEx.pl     blast AUGUSTUS protein sequences based on different extrinsic files and calculate
               a score for each prediction file to find optimal extrinsic parameter

SYNOPSIS

blastEx.pl [OPTIONS] --files=augustus1.gff,augustus2.gff

  --files                     list of AUGUSTUS prediction files in gff format with protein sequence  
  --best=output.out           contains file name of highest scoring extrinsic file
    
    
OPTIONS

    --help                          Print this help message
    --CPU                           Specifies the maximum number of CPUs that can be used during 
                                    computation
    --DB=path/to/blast/DB           Path to local blast database
    --dir=path/to/dir               Set path to working directory. In the working directory results
                                    and temporary files are stored
    --introns=hintsfile.gff         gff files with intron hints
    --log=lastEx.log                Log file
    --overwrite                     Overwrite existing files (except for species parameter files)

    

Format:
  seqname <TAB> source <TAB> feature <TAB> start <TAB> end <TAB> score <TAB> strand <TAB> frame <TAB> gene_id value <TAB> transcript_id value
                           

DESCRIPTION
      
  Example:

    blastEx.pl [OPTIONS] --files=augustus1.gff,augustus2.gff

ENDUSAGE



my $alpha = 1;              # factor for score (RNA-seq supported introns)
my $best_extrinsic;         # output file, contains name of best extrinsic file
my $beta = 1;               # factor for score (bases with protein support)
my $blastfactory;           # factory for BLAST
my $bool_good = "true";     # if true gene is good, if false, gene is bad
my $count_intron = 0;       # count the number of CDS per gene
my $CPU = 1;                # number of CPUs that can be used
my $DB;                     # path to local BLAST database
my $dir;                    # working superdirectory where programme is called from
my $false_pred = 0;         # number of falsely predicted introns (with no hint)
my @false_preds;            # array of all numbers of falsely predicted introns (for median)
my $good_mults = 0;         # all supported intron 'mult' entries summed up (all introns supported in gene)
my @gtf_files;              # gtf files from gff input files
my $help;                   # print usage   
my $highest_index;          # file index of the highest score
my $highest_score = '-inf'; # highest score 
my @ID;                     # gene ID
my $intron_mults = 0;       # all intron 'mult' entries summed up
my %introns;                # Contains information from intron file input. 
                            # Format $intron{seqname}{strand}[index]->{start} and ->{end}
my $introns;                # introns file name
my $log;                    # log file name
my $mults = 0;              # current 'mult' entries
my @no_blast;               # length of all sequences not blasted for each geneset (for median)
my $nr_of_bad = 0;          # number of bad genes
my $nr_of_genes = 0;        # number of all genes
my $nr_of_good = 0;         # number of good genes
my $nr_of_introns = 0;      # number of all introns in hints file
my $nr_of_introns_supp = 0; # number of all supported introns in prediction file
my $overwrite = 0;          # overwrite existing files (except for species parameter files)
my @parameters;             # parameters for blast
my @pred_files;             # list of input files in gff format
my $prog_dir;               # path to BLAST executable
my @prot_files;             # protein sequence files '.aa'
my %prot_seqs;              # hash of protein sequences
my $RNAseq_sup_introns = 0; # introns supported by RNA-seq data (sum of 'mult' entries)
my %scores;                 # score parts for each input file
my $score_bool = "true";    # true, if all scores are equal
my $standard = 0;           # index to standard value, if all scores are equal (currently extrinsic.m.0.1.b.1e1.cfg)
my $sup_mults = 0;          # sum of all supported 'mult' entries (all introns, but not necessarily all introns in gene are supported)
my $true_count = 0;         # counts the number of supported CDS per gene


if(@ARGV==0){
  print "$usage\n"; 
  exit(0);
}

GetOptions( 'best=s'      => \$best_extrinsic,
            'CPU=i'       => \$CPU,
            'DB=s'        => \$DB,
            'dir=s'       => \$dir,
            'files=s'     => \@pred_files,
            'introns=s'   => \$introns,
            'log=s'       => \$log,
            'overwrite!'  => \$overwrite,
            'param=f'     => \$alpha,
            'prog-dir=s'  => \$prog_dir,
            'help!'       => \$help);

if($help){
  print $usage;
  exit(0);
}

# if no working directory is set, use current directory
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

# parameters for blast without local database
@parameters = (-program    => "blastp",
               -db_name    => "nr",
               -prog_dir   => "$prog_dir",
               -expect     => '1e-5',
               -remote     => '1',
               -readmethod => "SearchIO");
if(defined($DB)){
  if(! -d $DB){
    print STDOUT "WARNING: Path to local BLAST database $DB does not exist. Programme will try to connect to server.\n";
    $blastfactory = Bio::Tools::Run::StandAloneBlastPlus->new(@parameters);
  }else{
    @parameters = (-program    => "blastp",
                   -db_name    => "nr",
                   -db_dir     => "$DB",
                   -prog_dir   => "$prog_dir",
                   -expect     => '1e-5',
                   -readmethod => "SearchIO");
    $blastfactory = Bio::Tools::Run::StandAloneBlastPlus->new(@parameters);
  }
}else{  
  $blastfactory = Bio::Tools::Run::StandAloneBlastPlus->new(@parameters);
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

    getAnnoFasta($pred_files[$i]);  # get protein sequence files if necessary  
    make_gtf($pred_files[$i]);      # maake gtf files if necessary
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
  $scores{$i}{"sup_introns"} = $sup_mults;  # introns supported by RNA-seq (log_10(mult))
  $scores{$i}{"good_genes"}  = $nr_of_good; # number of genes with full support (currently not part of the score)
  $scores{$i}{"bad_genes"}   = $nr_of_bad;  # number of genes without full support (currently not part of the score)
  $scores{$i}{"false_pred"}  = $false_pred; # falsely predicted introns in geneset i (no hint support)
  $scores{$i}{"false_neg"}   = $nr_of_introns - $nr_of_introns_supp; # false negative introns in geneset i (currently not part of the score)
  $scores{$i}{"pos"}         = 0;           # number of conserved (positives in stats section) in geneset i 
  $scores{$i}{"length"}      = 0;           # length of sequences blasted and in geneset i 
  $scores{$i}{"not_length"}  = 0;           # length of sequences blasted, but not in geneset i
  $scores{$i}{"not_blast"}   = 0;           # length of sequences in geneset i, but not blasted (with no blast hits)
  $scores{$i}{"not_pos"}     = 0;           # number of conserved (positives in stats section) not in geneset i  
  $scores{$i}{"nr_inhsp"}    = 0;           # number of high-scoring segment pairs in geneset i 
  $scores{$i}{"nr_outhsp"}   = 0;           # number of high-scoring segment pairs not in geneset i
}

# blast protein sequences
blast(); 
# print score hash to blastEx.stdout for manual checking     
print STDOUT Dumper(\%scores); 

my @score;
# calculation of median for falsely predicted introns and length of protein sequence without BLAST hit
foreach my $key (keys %scores){
  push(@false_preds, $scores{$key}{"false_pred"});
  push(@no_blast, $scores{$key}{"not_blast"});
}
my $median_fp = median(\@false_preds);
my $median_nb = median(\@no_blast);

# calculate score
for(my $i=0; $i<scalar(@gtf_files); $i++){
print LOG "\# ".(localtime).": supported introns ".$scores{$i}{"sup_introns"}." all intron mults $intron_mults positives ".$scores{$i}{"pos"}." length ".$scores{$i}{"length"}." length not in geneset ".$scores{$i}{"not_length"}." scores not blast ".$scores{$i}{"not_blast"}." good genes ".$scores{$i}{"good_genes"}." false neg ".$scores{$i}{"false_neg"}."\n";

  if($scores{$i}{"not_length"} == 0){
    $scores{$i}{"not_length"} = 1; 
  }
  if($scores{$i}{"length"} == 0){
    $scores{$i}{"length"} = 1;
  }
  my $score = $alpha * (($scores{$i}{"sup_introns"} - $scores{$i}{"false_pred"}) / $intron_mults) + ($scores{$i}{"sup_introns"} / ($scores{$i}{"sup_introns"} + $scores{$i}{"false_pred"})) + $alpha * ($scores{$i}{"length"} / ($scores{$i}{"not_length"} + $scores{$i}{"not_blast"})) + $beta *((($scores{$i}{"pos"} - $scores{$i}{"length"}) * $scores{$i}{"nr_outhsp"} / $scores{$i}{"length"}) - (($scores{$i}{"pos"} - $scores{$i}{"length"}) / ($scores{$i}{"not_pos"} - $scores{$i}{"not_length"})) - (($scores{$i}{"not_pos"} - $scores{$i}{"not_length"}) * $scores{$i}{"nr_inhsp"} / $scores{$i}{"not_length"})); 
  
  if($score >= $highest_score && $scores{$i}{"false_pred"} <= $median_fp && $scores{$i}{"not_blast"} <= $median_nb){
    $highest_score = $score;
    $highest_index = $i;
  }
  
  if($scores{$i}{"pos"} == 0 && $scores{$i}{"length"} == 0){
    print STDOUT "Number of positives and length is 0 for file $gtf_files[$i]. This is probably due to a BLAST connecting error. Please try again or use a local BLAST installation.\n";
  }
  push(@score, $score);
  print LOG "\# ".(localtime).": score $score, file $gtf_files[$i] and index $i\n";
  print STDOUT "\# ".(localtime).": score $score, file $gtf_files[$i] and index $i\n";
}
# check if all scores are equal
for(my $i=0; $i<scalar(@score)-1; $i++){
  if($score[$i] ne $score[$i+1]){
    $score_bool = "false";
  }
}
# if all scores are equal use standard parameters
if($score_bool eq "true"){
  $highest_index = $standard; # standard value, if all scores are equal
}

close(LOG) or die("Could not close file $log!\n");

# print best extrinsic file name to output
if(defined($best_extrinsic)){
  my @f = split(/\//, $gtf_files[$highest_index]);
  $f[-1] = "extrinsic".substr($f[-1],8,-3)."cfg"; #extrinsic.m.$malus[$i].b.$bonus[$j].cfg
  my $file = join("\/", @f);
  my $last_char = substr($file, -1);
  if($last_char eq "\/"){
    chop($file);
  }
  open (HIGHEST, ">$best_extrinsic") or die "Cannot open file: $best_extrinsic\n";
  print HIGHEST "$file";
  close (HIGHEST) or die("Could not close file $best_extrinsic!\n");
}

                           ############### sub functions ##############

# get protein sequences from files
sub get_prots{
  my $seqfile; 
  for(my $i=0; $i<scalar(@prot_files); $i++){
    print LOG "\# ".(localtime).": get protein sequences from $prot_files[$i]\n";
    $seqfile = Bio::SeqIO -> new (-file => "<$prot_files[$i]",
                                  -format => "fasta");
    while(my $seq = $seqfile -> next_seq){
      if(!defined($prot_seqs{$seq -> seq})){
        $prot_seqs{$seq -> seq}{"indices"} = "$i";
        $prot_seqs{$seq -> seq}{"sobj"} = $seq;
      }else{
        if($prot_seqs{$seq -> seq}{"indices"} !~ m/$i/){
          $prot_seqs{$seq -> seq}{"indices"} .= ",$i";
        }
      }
    }
  }
}

# blast sequences
sub blast{
  my $pm = new Parallel::ForkManager($CPU);
  foreach my $seq (keys %prot_seqs){
    my @file_index = split(/,/, $prot_seqs{$seq}{"indices"});
    $pm -> run_on_finish(
      sub{
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $scores_child) = @_;
        # update global hash
        foreach my $i (keys %scores){
          if(defined($scores_child->{$i}->{"pos"})){
            $scores{$i}{"pos"} += $scores_child->{$i}->{"pos"};
          }
          if(defined($scores_child->{$i}->{"length"})){
            $scores{$i}{"length"} += $scores_child->{$i}->{"length"};
          }
          if(defined($scores_child->{$i}->{"not_length"})){  
            $scores{$i}{"not_length"} += $scores_child->{$i}->{"not_length"};
          }
          if(defined($scores_child->{$i}->{"not_blast"})){
            $scores{$i}{"not_blast"} += $scores_child->{$i}->{"not_blast"};
          }
          if(defined($scores_child->{$i}->{"not_pos"})){
            $scores{$i}{"not_pos"} += $scores_child->{$i}->{"not_pos"};
          }
          if(defined($scores_child->{$i}->{"nr_inhsp"})){
            $scores{$i}{"nr_inhsp"} += $scores_child->{$i}->{"nr_inhsp"};
          }
          if(defined($scores_child->{$i}->{"nr_outhsp"})){
            $scores{$i}{"nr_outhsp"} += $scores_child->{$i}->{"nr_outhsp"};
          }
        }
      }
    );

    my $idx = 0;
    if(scalar(@file_index) != scalar(@prot_files) || scalar(@prot_files) == 1){
      $idx++;
      my $pid = $pm -> start($idx) and next;
      print STDOUT "\# ".(localtime).": blast ".$prot_seqs{$seq}{"sobj"} -> id."\n";
      my $result = $blastfactory -> run(-method => "blastp", -query => $prot_seqs{$seq}{"sobj"});
      print STDOUT "\# ".(localtime).": blast finished ".$prot_seqs{$seq}{"sobj"} -> id."\n";
      my %scores_c; my $size = keys(%{$result});
      # fill local hash
      if(%{$result}){
        print STDOUT "\# ".(localtime).": $idx get results for ".$result -> query_name()." number of hits ".$result -> num_hits."\nseq\n$seq\n";
        if($result -> num_hits != 0){
          my $hit = $result -> next_hit;
          if(defined($hit)){ 
            print STDOUT "\# ".(localtime).": hit for ".$result -> query_name()."\n";
            while(my $hsp = $hit -> next_hsp){ 
              print STDOUT "\# ".(localtime).": hsp for ".$result -> query_name()."\n";
              if(defined($prot_seqs{$seq})){
                my %hash; 
                foreach my $index (@file_index){
                  $hash{$index} = "";
                }
                for(my $i = 0; $i<scalar(@pred_files); $i++){
                  # current sequence is part of gene set $i
                  if(defined($hash{$i})){
                    print STDOUT "\# ".(localtime).": $i: num_conserved ".$hsp -> num_conserved." num_identical ".$hsp -> num_identical." length ".$prot_seqs{$seq}{"sobj"} -> length."\n";
                    $scores_c{$i}{"pos"} += $hsp -> num_conserved; 
                    $scores_c{$i}{"length"} += $prot_seqs{$seq}{"sobj"} -> length;
                    $scores_c{$i}{"nr_outhsp"}++;
                  # current sequence is not part of gene set $i
                  }else{
                    print STDOUT "\# ".(localtime).": $i: not length ".$prot_seqs{$seq}{"sobj"} -> length."\n";
                    $scores_c{$i}{"not_length"} += $prot_seqs{$seq}{"sobj"} -> length;
                    $scores_c{$i}{"not_pos"} += $hsp -> num_conserved;
                    $scores_c{$i}{"nr_inhsp"}++;
                  }
                }
              }
            }
          }
        }else{
          for(my $j = 0; $j<scalar(@file_index); $j++){
            $scores_c{$file_index[$j]}{"not_blast"} += $prot_seqs{$seq}{"sobj"} -> length;
            $scores_c{$file_index[$j]}{"nr_outhsp"}++;
          }
          print STDOUT "\# ".(localtime).": no hits ".$prot_seqs{$seq}{"sobj"} -> id."\n";
        }
      }else{
        for(my $j = 0; $j<scalar(@file_index); $j++){
          $scores_c{$file_index[$j]}{"not_blast"} += $prot_seqs{$seq}{"sobj"} -> length;
          $scores_c{$file_index[$j]}{"nr_outhsp"}++;
        }
        print STDOUT "\# ".(localtime).": no result ".$prot_seqs{$seq}{"sobj"} -> id."\n";
      }
      $pm -> finish(0, \%scores_c);
    }else{
      print STDOUT "\# ".(localtime).": not blast ".$prot_seqs{$seq}{"sobj"} -> id."\n";
    }
  }
  $pm -> wait_all_children;
}

# read in introns
sub introns{
  open (INTRONS, $introns) or die "Cannot open file: $introns\n";
  print LOG "\# ".(localtime).": read in introns from $introns\n";
  while(<INTRONS>){
    chomp;
    my @line = split(/\t/, $_);
    if(scalar(@line) == 9){
      $introns{$line[0]}{$line[6]}{$line[3]}{$line[4]} = log10(1 + $line[5]);
      $nr_of_introns++;
      $intron_mults += log10(1 + $line[5]);
    }
  }  
  close(INTRONS) or die("Could not close file $introns!\n");
}

# call getAnnoFasta.pl
sub getAnnoFasta{
  my $AUG_pred = shift;
  @_ = split(/\//, $AUG_pred);
  my $name_base = substr($_[-1],0,-4);
  my $prot_file = substr($AUG_pred,0,-4).".aa";
  if(!uptodate([$AUG_pred],[$prot_file])  || $overwrite){
    my $string = find("getAnnoFasta.pl");
    my $errorfile = "$dir/errors/getAnnoFasta.$name_base.stderr";
    my $perlCmdString = "perl $string $AUG_pred 2>$errorfile";
    print LOG "\# ".(localtime).": make a fasta file with protein sequences for $AUG_pred\n";
    print LOG "$perlCmdString\n\n";
    system("$perlCmdString")==0 or die("failed to execute: $!\n");
  }
  push(@prot_files, $prot_file); 
}

# make gtf file
sub make_gtf{
  my $AUG_pred = shift;
  @_ = split(/\//, $AUG_pred);
  my $name_base = substr($_[-1],0,-4);
  my $gtf_file = substr($AUG_pred,0,-4).".gtf";
  if(!uptodate([$AUG_pred],[$gtf_file])  || $overwrite){
    my $errorfile = "$dir/errors/gtf2gff.$name_base.stderr";
    my $perlstring = find("gtf2gff.pl");
    my $cmdString = "cat $AUG_pred | perl -ne 'if(m/\\tAUGUSTUS\\t/){print \$_;}' | perl $perlstring --printExon --out=$gtf_file 2>$errorfile";
    print "$cmdString\n\n";
    print LOG "\# ".(localtime).": make a gtf file from $AUG_pred\n";
    print LOG "$cmdString\n\n";
    system("$cmdString")==0 or die("failed to execute: $!\n");
  }
  push(@gtf_files, $gtf_file); 
}

# compare input files to introns from hints file
sub compare_introns{
  my $gtf_file = shift;
  my $prev_ID = "no_ID";

  $sup_mults = 0;
  $good_mults = 0;
  $false_pred = 0;
  $nr_of_genes = 0;
  $count_intron = 0;
  $nr_of_good = 0;
  $nr_of_introns_supp = 0;
  $nr_of_bad = 0;

  open (GTF, "<".$gtf_file) or die "Cannot open file: $gtf_file\n";
  while(<GTF>){
    chomp;
    my @line = split(/\t/, $_);
    if(scalar(@line) == 9){
      if(($line[2] ne "start_codon") && ($line[2] ne "stop_codon") && ($line[2] ne "intron")){
        next;
      }else{
        @ID = split(/\s/,$line[8]); # gene_id "gene:SPAC212.11"; transcript_id "transcript:SPAC212.11.1";
        # new gene starts
        if($prev_ID ne $ID[1]){
          if($count_intron != 0){
            $nr_of_genes++;
            reset_param();
          }
        }
        if($line[2] eq "intron"){
          # check if intron is in hash made of intron input
          if(defined($introns{$line[0]}{$line[6]}{$line[3]}{$line[4]})){
            $true_count++;
            $nr_of_introns_supp++;
            $mults += $introns{$line[0]}{$line[6]}{$line[3]}{$line[4]};
          }else{
            $false_pred++;
          }

          $count_intron++;
        }
        $prev_ID = $ID[1];
      }
    }
  }
  $nr_of_genes++;
  reset_param(); 

  close (GTF) or die("Could not close file $gtf_file!\n");
}

# reset parameters
sub reset_param{
  if($true_count != $count_intron && $count_intron != 0){
    $bool_good = "false";
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
  $count_intron = 0;
  $true_count = 0;
  $bool_good = "true";
  $mults = 0;
}
