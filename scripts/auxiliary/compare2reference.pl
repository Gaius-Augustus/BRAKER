#!/usr/bin/perl

####################################################################################################
#                                                                                                  #
# compare2reference.pl - compares prediction file to reference and get some statistics about       #
#                        predicted genes, like how many genes are correctly predicted, have        #
#                        wrong start or stop codons, are wrongly splitted or connected or          #
#                        were not predicted at all                                                 #
#                        can also create link files for each case                                  #
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
# | first outline                   | Simone Lange   |22.04.2015 |
# | add link file outputs           |                |           |
# ----------------------------------------------------------------
 
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(ceil);
use Data::Dumper;



my $usage = <<'ENDUSAGE';

compare2reference.pl     compare prediction file to reference and get statistics (input files need 'gene' lines

SYNOPSIS

compare2reference.pl [OPTIONS] --in=prediction.gff3 --ref=reference.gff3

  prediction.gff3         prediction file in gff, gtf or gff3 format (needs 'gene' line entries)
  reference.gff3          reference file in gff, gtf or gff3 format (needs 'gene' line entries)
  
    
OPTIONS

    --help                          Print this help message
    --in=prediction.gff3            Prediction file (needs 'mRNA' line entries)
    --JBROWSE_PATH=/path/           Set path to JBrowse (if not specified as environment variable).
    --ref=reference.gff3            Reference file (needs 'gene' line entries)
    --tracks=track1,track2          List of track names. Separators are ',' or spaces


Format:
  seqname <TAB> source <TAB> feature <TAB> start <TAB> end <TAB> score <TAB> strand <TAB> gene_id value <TAB> link
                           

DESCRIPTION
      
  Example:

    compare2reference.pl [OPTIONS] --in=prediction.gff3 --ref=reference.gff3

ENDUSAGE


my $all_genes_in = 0;           # count number of all genes in prediction file
my $all_genes_ref = 0;          # count number of all genes in reference file
my $bool_match = "false";       # true if transcript has right start and stop codon and supported exons, false if only start and stop codon is right, check if CDS lines of current mRNA should be checked for matching
my $bool_missing = "false";     # true if transcript is missing, i.e. transcript in reference but not in prediction
my $bool_not_split = "no_case"; # true if transcript is not falsely splitted, false if transcript is falsely splitted, no case if this case should not be considered
my $bool_wrong_con = "false";   # true if transcript is falsely connected
my $bool_wrong_start = "false"; # true if transcript has a wrong start codon
my $bool_wrong_stop = "false";  # true if transcript has a wrong stop codon
my %CDS;                        # hash of arrays of prediction input (format $CDS{seqname}{strand}{start}{end} = ID
my $cds = 0;                    # counts number of transcripts with matching start and stop codon but at least one wrong CDS entry
my %CDS_count;                  # counts the number of CDS entries for each transcript
my $CDS_trans = 0;              # number of CDS in transcript (reference)
my %current_IDS;                # IDs of possible matching genes in prediction with counter
my $help;                       # print usage 
my $in;                         # prediction file name
my %in_ends;                    # hash of arrays of prediction ends (format: $in_end{seqname}{strand}{end} = start)
my %in_starts;                  # hash of arrays of prediction starts (format: $in_start{seqname}{strand}{start} = end)
my $jbrowse_path;               # JBrowse path, higher priority than $PATH_TO_JBROWSE
my $link_start;                 # link base start
my $link_end;                   # link base end
my $match = 0;                  # counts number of matching genes
my $missing_gene = 0;           # counts number of missing genes (i.e. not predicted, substract wrongly predicted genes at the end)
my $mRNA_check_join;            # important columns of mRNA line that is currently checked for matching
my $mRNA_check_link;            # link for mRNA line that is currently checked for matching
my $PATH_TO_JBROWSE = $ENV{'JBROWSE_PATH'}; # export JBROWSE_PATH=/var/www/jbrowse/JBrowse-1.11.5/ 
my @tracks;                     # list of track labels
my $ref;                        # reference file name
my %ref_end;                    # hash of arrays of reference input (format $ref{seqname}{strand}{end} = start
my %ref_start;                  # hash of arrays of reference input (format $ref{seqname}{strand}{start} = end
my $wrong_connect = 0;          # count number of wrongly connected genes
my $wrong_start = 0;            # count number of genes with wrong start codon (without wrongly splitted/connected genes)
my $wrong_stop = 0;             # count number of genes with wrong start codon (without wrongly splitted/connected genes)
my $wrong_split = 0;            # count number of wrongly splitted genes  
my $wrong_pred = 0;             # count number of wrongly predicted genes
my $out_start;                  # output file name for transcripts with wrong start codon
my $out_start_stop;             # output file name for transcripts with matching start and stop codons but at least one wrong CDS entry
my $out_stop;                   # output file name for transcripts with wrong stop codon
my $out_connect;                # output file name for falsely connected transcripts 
my $out_split;                  # output file name for falsely splitted transcripts 
my $out_match;                  # output file name for matching transcripts 
my $out_miss;                   # output file name for missing or falsely predicted (from reference point of view) transcripts 
my $out_wpred;                  # output file name for falsely predicted (from prediction point of view) transcripts transcripts 


if(@ARGV==0){
  print "$usage\n"; 
  exit(0);
}

GetOptions( 'in=s'            => \$in,
            'JBROWSE_PATH=s'  => \$jbrowse_path,
            'ref=s'           => \$ref,
            'tracks=s'        => \@tracks, 
            'help!'           => \$help);

if($help){
  print $usage;
  exit(0);
}

if(defined($jbrowse_path)){
  $PATH_TO_JBROWSE = $jbrowse_path;
}

if(!$ENV{'JBROWSE_PATH'} && !defined($jbrowse_path)){
  print STDOUT "WARNING: No path for JBrowse was set. Programme will use '/var/www/jbrowse/JBrowse-1.11.5'. If you want to use a different path, please export an environment variable or use --JBROWSE_PATH=path/to/jbrowse.\n";
  $PATH_TO_JBROWSE = "/var/www/jbrowse/JBrowse-1.11.5/";
}

if($PATH_TO_JBROWSE !~ /\/$/){
  $PATH_TO_JBROWSE .= "/";
}

if($PATH_TO_JBROWSE !~ /\/^/){
  $PATH_TO_JBROWSE = "/$PATH_TO_JBROWSE";
}

$link_start = "file://".$PATH_TO_JBROWSE."index.html?loc=";

if(@tracks){
  @tracks = split(/[\s,]/, join(',',@tracks));
  my $join = join("%2C",@tracks);
  $link_end = "&tracks=".$join."&highlight=";
}else{
  $link_end = "";
}

# set $ref
if(!defined($ref)){
  $ref = $ARGV[0];
}

# set $ins
if(!defined($in)){
  $in = $ARGV[1];
}

# check whether the input files exists
if(! -f $in || ! -f $ref){
    print STDERR "Needs both reference and input file, but at least one of $in and $ref does not exist. Please check.\n";
    exit(1);
}else{
  $out_start = substr($in,0,-4)."wrong.start.html";
  $out_start_stop = substr($in,0,-4)."wrong.CDS.html";
  $out_stop = substr($in,0,-4)."wrong.stop.html";
  $out_connect = substr($in,0,-4)."wrong.connect.html";
  $out_split = substr($in,0,-4)."wrong.split.html";
  $out_match = substr($in,0,-4)."match.html";
  $out_miss = substr($in,0,-4)."miss.html";
  $out_wpred = substr($in,0,-4)."wrong_miss.pred.html";

  open(START, ">", $out_start) or die("Could not open file $out_start!\n");
  print START "<html>\n";
  print START "<table>\n";
  open(CDS, ">", $out_start_stop) or die("Could not open file $out_start_stop!\n");
  print CDS "<html>\n";
  print CDS "<table>\n";
  open(STOP, ">", $out_stop) or die("Could not open file $out_stop!\n");
  print STOP "<html>\n";
  print STOP "<table>\n";
  open(CON, ">", $out_connect) or die("Could not open file $out_connect!\n");
  print CON "<html>\n";
  print CON "<table>\n";
  open(SPLIT, ">", $out_split) or die("Could not open file $out_split!\n");
  print SPLIT "<html>\n";
  print SPLIT "<table>\n";
  open(MATCH, ">", $out_match) or die("Could not open file $out_match!\n");
  print MATCH "<html>\n";
  print MATCH "<table>\n";
  open(MISS, ">", $out_miss) or die("Could not open file $out_miss!\n");
  print MISS "<html>\n";
  print MISS "<table>\n";
  open(PRED, ">", $out_wpred) or die("Could not open file $out_wpred!\n");
  print PRED "<html>\n";
  print PRED "<table>\n";
  get_ref();
  input();  
  compare();
  print START "<Xtable>\n";
  print START "</html>\n";
  print CDS "<Xtable>\n";
  print CDS "</html>\n";
  print STOP "<Xtable>\n";
  print STOP "</html>\n";
  print CON "<Xtable>\n";
  print CON "</html>\n";
  print SPLIT "<Xtable>\n";
  print SPLIT "</html>\n";
  print MATCH "<Xtable>\n";
  print MATCH "</html>\n";
  print MISS "<Xtable>\n";
  print MISS "</html>\n";
  print PRED "<Xtable>\n";
  print PRED "</html>\n";
  close (START) or die("Could not close file $out_start!\n");
  close (CDS) or die("Could not close file $out_start_stop!\n");
  close (STOP) or die("Could not close file $out_stop!\n");
  close (CON) or die("Could not close file $out_connect!\n");
  close (SPLIT) or die("Could not close file $out_split!\n");
  close (MATCH) or die("Could not close file $out_match!\n");
  close (MISS) or die("Could not close file $out_miss!\n");
  close (PRED) or die("Could not close file $out_wpred!\n");  
    
}

print "Number of wrong start codons:\t$wrong_start\t".(($wrong_start / $all_genes_in) * 100)."%\n";
print "Number of wrong stop codons:\t$wrong_stop\t".(($wrong_stop / $all_genes_in) * 100)."%\n";
print "Number of wrongly splitted transcripts:\t$wrong_split\t".(($wrong_split / $all_genes_in) * 100)."%\n";
print "Number of wrongly connected transcripts:\t$wrong_connect\t".(($wrong_connect / $all_genes_in) * 100)."%\n";
print "Number of wrongly predicted transcripts:\t$wrong_pred\t".(($wrong_pred / $all_genes_in) * 100)."%\n";
print "Number of wrongly predicted and/or missing transcripts:\t".($missing_gene)."\t".(($missing_gene / $all_genes_in) * 100)."%\n";
print "Number of transcripts with matching start and stop codon but at least one wrong CDS entry:\t$cds\t".(($cds / $all_genes_in) * 100)."%\n";
print "Number of matching transcripts:\t$match\t".(($match / $all_genes_in) * 100)."%\n";
print "Number of all transcripts (prediction):\t$all_genes_in\t".($all_genes_in / $all_genes_ref)."\n";
print "Number of all transcripts (reference):\t$all_genes_ref\n";


                           ############### sub functions ##############
# get the reference genes for cross checking
sub get_ref{
  open (REF, "<".$ref) or die "Cannot open file: $ref\n";
  while(<REF>){
    chomp;
    my @line = split(/\t/, $_);
    if(scalar(@line) == 9){
      if($line[2] eq "mRNA"){
        push(@{$ref_end{$line[0]}{$line[6]}{$line[4]}},  $line[3]) unless grep{$_ == $line[3]} @{$ref_end{$line[0]}{$line[6]}{$line[4]}};
        push(@{$ref_start{$line[0]}{$line[6]}{$line[3]}},  $line[4]) unless grep{$_ == $line[4]} @{$ref_start{$line[0]}{$line[6]}{$line[3]}};
      }
    }
  }  
  close (REF) or die("Could not close file $ref!\n");
}

# read in prediction file
sub input{
  open (IN, $in) or die "Cannot open file: $in\n";
  while(<IN>){
    chomp;
    my @line = split(/\t/, $_);
    if(scalar(@line) == 9){
      if($line[2] eq "mRNA"){
        $all_genes_in++;
        push(@{$in_ends{$line[0]}{$line[6]}{$line[4]}},  $line[3]) unless grep{$_ == $line[3]} @{$in_ends{$line[0]}{$line[6]}{$line[4]}};
        push(@{$in_starts{$line[0]}{$line[6]}{$line[3]}},  $line[4]) unless grep{$_ == $line[4]} @{$in_starts{$line[0]}{$line[6]}{$line[3]}};
        if(!defined($ref_start{$line[0]}{$line[6]}{$line[3]}) && !defined($ref_end{$line[0]}{$line[6]}{$line[4]})){
          my $link_pos = "$line[0]%3A$line[3]..$line[4]";
          my $link = $link_start.$link_pos.$link_end; 
          print PRED "<tr>\n";
          splice(@line,7,1);
          my $join = join("<td>",@line);
          print PRED "<td>".$join."<td><a href = $link>link</a> </td>\n";
          $wrong_pred++;
        } 
      }
      # store CDS entries for checking if gene matches not only start and stop codon
      if($line[2] eq "CDS"){
        @_ = split(/;/, $line[8]);
        push(@{$CDS{$line[0]}{$line[6]}{$line[3]}{$line[4]}},  $_[1]) unless grep{$_ eq $_[1]} @{$CDS{$line[0]}{$line[6]}{$line[3]}{$line[4]}};
        if(defined($CDS_count{$_[1]})){
          $CDS_count{$_[1]}++;
        }else{
          $CDS_count{$_[1]} = 1;
        }
      }
    }
  }  
  close (IN) or die("Could not close file $in!\n");
}



# compare reference file to prediction file
sub compare{
  open (REF, "<".$ref) or die "Cannot open file: $ref\n";
  while(<REF>){
    chomp;
    my @line = split(/\t/, $_);
    if(scalar(@line) == 9){
      if($line[2] eq "mRNA"){
        # check for last mRNA line
        if($bool_match eq "check"){
          # check for each possible matching transcript if all CDS entries are right and if the number of all CDS entries of those transcript are also matching 
          foreach my $ID (keys %current_IDS){
            if($current_IDS{$ID} == $CDS_trans && $CDS_count{$ID} == $CDS_trans){
              $bool_match = "true";
              last;
            }
          }
          if($bool_match eq "true"){
            print MATCH "<tr>\n";
            print MATCH "<td>".$mRNA_check_join."<td><a href = $mRNA_check_link>link</a> </td>\n";
            $match++;
          }else{
            print CDS "<tr>\n";
            print CDS "<td>".$mRNA_check_join."<td><a href = $mRNA_check_link>link</a> </td>\n";
            $cds++;
          }
          $bool_match = "false";
          %current_IDS = ();
          $CDS_trans = 0;
        }
        # check for new mRNA line
        $all_genes_ref++;
        my $link_pos = "$line[0]%3A$line[3]..$line[4]";
        splice(@line,7,1);
        my $join = join("<td>",@line);
        my $link = $link_start.$link_pos.$link_end; 
        # matching stop codon
        if(defined($in_ends{$line[0]}{$line[6]}{$line[4]}) && !defined($in_starts{$line[0]}{$line[6]}{$line[3]})){
          foreach my $start (@{$in_ends{$line[0]}{$line[6]}{$line[4]}}){
            if(defined($ref_start{$line[0]}{$line[6]}{$start})){
              $bool_wrong_con = "true";
              foreach my $end (@{$ref_start{$line[0]}{$line[6]}{$start}}){
                if($end == $line[4]){
                  $bool_missing = "true";
                  $bool_wrong_con = "false";
                  last;
                }
              }
              last;
            }else{
              $bool_wrong_start = "true";
            }
          }
        # matching start codon
        }elsif(!defined($in_ends{$line[0]}{$line[6]}{$line[4]}) && defined($in_starts{$line[0]}{$line[6]}{$line[3]})){
          foreach my $end (@{$in_starts{$line[0]}{$line[6]}{$line[3]}}){
            if(defined($ref_end{$line[0]}{$line[6]}{$end})){
              $bool_wrong_con = "true";            
              foreach my $start (@{$ref_end{$line[0]}{$line[6]}{$end}}){
                if($start == $line[3]){
                  $bool_missing = "true";
                  $bool_wrong_con = "false";
                  last;
                }
              }
              last;
            }else{
              $bool_wrong_stop = "true";
            }
          }
        # matching start and stop codon
        }elsif(defined($in_ends{$line[0]}{$line[6]}{$line[4]}) && defined($in_starts{$line[0]}{$line[6]}{$line[3]})){
          $bool_not_split = "false";
          foreach my $start (@{$in_ends{$line[0]}{$line[6]}{$line[4]}}){
            if($start == $line[3]){
              $bool_not_split = "true";
              last;
            }
          }
          foreach my $end (@{$in_starts{$line[0]}{$line[6]}{$line[3]}}){
            if($end == $line[4]){
              $bool_not_split = "true";
              last;
            }
          }
          if($bool_not_split eq "true"){
            $bool_match = "check"; 
            $mRNA_check_join = $join;
            $mRNA_check_link = $link;  
          }
        # neither start nor stop codon matching
        }elsif(!defined($in_ends{$line[0]}{$line[6]}{$line[4]}) && !defined($in_starts{$line[0]}{$line[6]}{$line[3]})){
          $bool_missing = "true";
        }

        if($bool_wrong_start eq "true"){
          print START "<tr>\n";
          print START "<td>".$join."<td><a href = $link>link</a> </td>\n";
          $wrong_start++;
        }elsif($bool_wrong_stop eq "true"){
          print STOP "<tr>\n";
          print STOP "<td>".$join."<td><a href = $link>link</a> </td>\n";
          $wrong_stop++;
        }elsif($bool_wrong_con eq "true"){
          print CON "<tr>\n";
          print CON "<td>".$join."<td><a href = $link>link</a> </td>\n";
          $wrong_connect++;
        }
        if($bool_not_split eq "false"){
          print SPLIT "<tr>\n";
          print SPLIT "<td>".$join."<td><a href = $link>link</a> </td>\n";
          $wrong_split++;
        }
        if($bool_missing eq "true"){
          print MISS "<tr>\n";
          print MISS "<td>".$join."<td><a href = $link>link</a> </td>\n";
          $missing_gene++;
        }
        $bool_wrong_con = "false";
        $bool_missing = "false";
        $bool_wrong_start = "false";
        $bool_wrong_stop = "false";
        $bool_not_split = "no_case";
      }elsif($line[2] eq "CDS" && $bool_match eq "check"){     
        foreach my $ID (@{$CDS{$line[0]}{$line[6]}{$line[3]}{$line[4]}}){
          if($CDS_trans == 0){
            $current_IDS{$ID} = 1;
          }else{    
            if(defined($current_IDS{$ID})){
              $current_IDS{$ID}++;
            }
          }
        }
        $CDS_trans++;
      }
    }
  }
  close (REF) or die("Could not close file $ref!\n");
}


