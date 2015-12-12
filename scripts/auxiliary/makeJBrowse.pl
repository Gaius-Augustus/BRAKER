#!/usr/bin/perl
####################################################################################################
#                                                                                                  #
# makeJBrowse.pl - make html files with links to JBrowse on given coordinates                      #
#                  each file contains only links to positions that two files have in common        #
#                  e.g. makeBrowser.0.1.html contains all genes that only these two files          #
#                  have in common. 'filenames.map' contains the file names and their number        #
#                  in tabulator separated format                                                   #
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
# | creation of file                | Simone Lange   |09.12.2014 |
# ----------------------------------------------------------------
# 
use strict;
use warnings;
use Getopt::Long;


my $usage = <<'ENDUSAGE';

makeJBrowse.pl     make html files with links to JBrowse on given coordinates

SYNOPSIS

makeJBrowse.pl [OPTIONS] --in=file1.gff3,file2.gff3

  --in=file1.gff3,file2.gff3          list of files in gff, gtf or gff3 format. Separators are ',' or spaces
  
    
OPTIONS

    --help                          Output this help message
    --tracks=track1,track2          List of track names. Separators are ',' or spaces
    --JBROWSE_PATH=/path/           Set path to JBrowse (if not specified as environment variable).


Format:
  seqname <TAB> source <TAB> feature <TAB> start <TAB> end <TAB> score <TAB> strand <TAB> frame <TAB> gene_id value <TAB> transcript_id value <TAB> link
                           

DESCRIPTION
      
  Example:

    makeJBrowse.pl [OPTIONS] --in=ref.gff3,augustus.gff3

ENDUSAGE


my @hashes;           # array of hash references for each input file: $hash{seqname}{strand}{start}{end}
my $help;             # print usage
my @in;               # list of input files in gff, gtf or gff3 format. Separators are ',' or spaces
my $jbrowse_path;     # JBrowse path, higher priority than $PATH_TO_JBROWSE
my $link_start;       # link base start
my $link_end;         # link base end
my @out_files;        # output file names array
my %out_files;        # output file names hash
my $PATH_TO_JBROWSE = $ENV{'JBROWSE_PATH'}; # export JBROWSE_PATH=/var/www/jbrowse/JBrowse-1.11.5/
my @tracks;           # list of track labels

if(@ARGV==0){
  print "$usage\n"; 
  exit(0);
}

GetOptions( 'in=s'            => \@in,
            'JBROWSE_PATH=s'  => \$jbrowse_path,
            'tracks=s'        => \@tracks,
            'help!'           => \$help);

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
  print STDOUT "WARNING: No track labels were assigned.\n" 
}


if($help){
  print $usage;
  exit(0);
}

if(@in){
  @in = split(/[\s,]/, join(',',@in));
  for(my $i=0; $i<scalar(@in); $i++){
    if(! -e $in[$i]){
      print STDERR "ERROR: Input file $in[$i] does not exist. Please check.\n";
      exit(1);
    }
  }
  if(scalar(@in) < 2){
    print STDERR "ERROR: You have to assign at least two files. Please check.\n";
    exit(1);
  }else{
    readIn();
    compareFiles();
  }
}


                           ############### sub functions ##############

# read in files and store them in hash
sub readIn{
  open (MAP, ">", "filenames.map") or die "Cannot open file: filenames.map\n";
  for(my $i=0; $i<scalar(@in); $i++){
    open (IN, $in[$i]) or die "Cannot open file: $in[$i]\n";
    print MAP "$i\t$in[$i]\n";
    my %hash = ();
    while(<IN>){
      chomp;
      my @line = split(/\t/, $_);
      $hash{$line[0]}{$line[6]}{$line[3]}{$line[4]} = "";
    }
    push(@hashes, \%hash);
    close IN;
  }
  close MAP;
}


# compare the input files and only save entries in output files that only two of them have in common
sub compareFiles{
  for(my $i=0; $i<scalar(@in)-1; $i++){
    open (IN, $in[$i]) or die "Cannot open file: $in[$i]\n";
    while(<IN>){
      chomp;
      my @line = split(/\t/, $_);
      my $def_count = 0;
      next if ($line[2] eq "start_codon" || $line[2] eq "stop_codon");
      for(my $j=0; $j<scalar(@hashes); $j++){ 
        next if ($i == $j);
        if(defined($hashes[$j]{$line[0]}{$line[6]}{$line[3]}{$line[4]})){
          $def_count++;
          $hashes[$i]{$line[0]}{$line[6]}{$line[3]}{$line[4]} = $j;
        }
      }
      if(($def_count == 1 && $i < $hashes[$i]{$line[0]}{$line[6]}{$line[3]}{$line[4]})){
        my $link_pos = "$line[0]%3A$line[3]..$line[4]";
        my $link = $link_start.$link_pos.$link_end; 
        my $out = "makeBrowser.$i.".$hashes[$i]{$line[0]}{$line[6]}{$line[3]}{$line[4]}.".html";
        if(!defined($out_files{$out})){
          $out_files{$out} = "";
          open(OUT, ">", $out) or die("Could not open file $out!\n");
          print OUT "<html>\n";
          print OUT "<table>\n";
          push(@out_files, $out);         
        }else{
          open(OUT, ">>", $out) or die("Could not open file $out!\n");
        }
        print OUT "<tr>\n";
        my $join = join("<td>",@line);
        print OUT "<td>".$join."<td><a href = $link>link</a> </td>\n";
        close(OUT);
      }
    }
    close(IN);
  }
  foreach my $out (@out_files){
    open(OUT, ">>", $out) or die("Could not open file $out!\n");  
    print OUT "<Xtable>\n";
    print OUT "</html>\n";
  }
}


