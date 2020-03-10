#!/usr/bin/env perl

####################################################################################################
#                                                                                                  #
# merge_transcript_sets.pl   - merge N transcript sets in such a way that transcripts that were    #
#                              missing in the files up to N (1..N-1) are added to the set of       #
#                              transcripts.                                                        #
#                              Attention: nonredundant transcript names over all sets are assumed! #
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
# | file creation                     | Katharina Hoff |18.12.2019 |
# ------------------------------------------------------------------

use strict;
use warnings;
use Getopt::Long;

my $usage = <<'ENDUSAGE';

merge_transcript_sets.pl  merge N transcript sets in such a way that transcripts that were
                          missing in the files up to N (1..N-1) are added to the set of
                          transcripts.
                          Attention: nonredundant transcript names over all sets are assumed!

SYNOPSIS

merge_transcript_sets.pl [OPTIONS] set1.gtf set2.gtf ...

  setN.gtf               file with gene predictions in gtf format


OPTIONS

    --help               print this help message


DESCRIPTION

  Example:

    merge_transcript_sets.pl set1.gtf set2.gtf > out.gtf

ENDUSAGE

my $help;

GetOptions( 'help!'     => \$help);

if( $help || (scalar(@ARGV) < 1) ) {
  print $usage;
  exit(0);
}

my %txid_to_elements;
my %uniq_struct_to_txid;
my $file_counter = 1;
foreach(@ARGV){
    my $file = $_;
    print STDERR "Processing file $file\n";
    # sort gtf file (just to be sure)
    my $cmdString = "cat $file | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -n -k 3,3 | sort -s -k 1,1 > $file".".sorted";
    system($cmdString) == 0 or die("ERROR in file " . __FILE__ ." at line "
        . __LINE__ ."\nFailed to execute: $cmdString!\n");
    open(IN, "<", $file.".sorted") or die ("ERROR: in file " . __FILE__ ." at line "
        . __LINE__ ."\nFailed to open file $file".".sorted for reading!\n");
    my @store_for_txid = (); # gene lines do not contain transcript_id and can only be appended later
    my %txid_to_struct_local;
    while(<IN>){
        if(not($_=~m/\#/)){
            my $line = $_;
            $line =~ s/\"([^"]+)\"/\"file_${file_counter}_$1\"/g;
            my $txid;
            if($line =~ m/transcript_id/){
                $line =~ m/transcript_id "([^"]+)";/;
                $txid = $1;
                push(@{$txid_to_elements{$txid}}, $line);
                foreach(@store_for_txid){
                    push(@{$txid_to_elements{$txid}}, $_)
                }
                @store_for_txid = ();
            }else{
                $line =~ s/\t([\t]+)$/\tfile_${file_counter}_$1/;
                push(@store_for_txid, $line);
            }
            # currently, UTR features are ignored
            if($line =~ m/\tCDS\t/){
                my @t = split(/\t/, $line);
                if(not(defined($txid_to_struct_local{$txid}))){
                    $txid_to_struct_local{$txid} = $t[0]."_".$t[3]."_".$t[4]."_".$t[6];
                }else{
                    $txid_to_struct_local{$txid} .= "_".$t[0]."_".$t[3]."_".$t[4]."_".$t[6];
                }
            }
        }
    }
    close(IN) or die ("ERROR: in file " . __FILE__ ." at line "
        . __LINE__ ."\nFailed to close file $file".".sorted!\n");
    unlink($file.".sorted") or die ("ERROR: in file " . __FILE__ ." at line "
        . __LINE__ ."\nFailed to delete file $file".".sorted!\n");
    # always keep the first occuring transcript structure, only add from other gene sets if it has not been in the set, yet
    # this might discard alternative UTR splicing isoforms at present
    while (my ($key, $value) = each (%txid_to_struct_local)){
        #print "key is $key and value is $value\n";
        if(not(defined($uniq_struct_to_txid{$value}))){
            #print "adding transcript\n";
            $uniq_struct_to_txid{$value} = $key;
        }
    }
    $file_counter++;
}

# print result
while (my ($key, $value) = each (%uniq_struct_to_txid)){
    #print "\n".$key."\n";
    foreach(@{$txid_to_elements{$value}}){
        print $_;
    }    
}