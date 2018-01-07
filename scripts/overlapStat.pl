#!/usr/bin/perl
# overlapStat.pl
# Examine the combined overlap of a set of lists
# or a set of integer intervals like different genome annotations.
#
# Mario Stanke, 7.6.2007
#

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $usage = "Usage:\n";
$usage .= "   overlapStat.pl set1 set2 set3 ... setn\n";
$usage .= "   Examine the combined overlap of a set of n sets.\n";
$usage .= "   Computes the sizes of all 2^n-1 logical combinations like, e.g.\n";
$usage .= "   'in set 1 and in set 2 and not in set 3'.\n";
$usage .= "Options:\n";
$usage .= "      --outfiles    If this flag is set, combination sets are created.\n";
$usage .= "   The n input files can have one of two formats.\n";
$usage .= "   Firstly, they can be a simple list of names:\n";
$usage .= "   set1:                  set2:             ...\n";
$usage .= "      Daniel                 Gerlinde\n";
$usage .= "      Lars                   Irfan\n";
$usage .= "      Stephanie              Stephanie\n";
$usage .= "      Noah                   Lars\n";
$usage .= "      Nikolai             \n";
$usage .= "      Mario               \n";
$usage .= "   Secondly, they can be labeled integer intervals:\n";
$usage .= "   set1:                  set2:                    ...\n";
$usage .= "      February  3  12        January   1  23\n";
$usage .= "      March     15 17        February  5  10\n";
$usage .= "      April     5  22        April     20 25\n";
$usage .= "      April     26 27        \n";
$usage .= "   The second format can in particular be used to compare sets of genome annotations.\n";
$usage .= "   For example, as a special case you can compute base level sensitivity and specificity.\n";
$usage .= "   The input sets do not need to be sorted. Repetitions of the same name and\n";
$usage .= "   self-overlapping of a set are allowed.\n";
my $outfiles = 0; # if true, output the files with the list in addition to the stats
GetOptions('outfiles!'=>\$outfiles);

my $n = scalar @ARGV;
if ($n < 1) {
    print $usage;
    exit;
}
my $simple = 0;
my @rawsets;    # raw: as given
my @sets;       # broken in pieces given by the canonical breakpoints
my %alllabels = ();
my %union = (); # contains the interval lengths for items

#
# data structure:
# rawsets is array of 
# hashes with the labels are the keys
# and where the values are references 
# to lists of references to arrays
# e.g. $rawsets[0]->{"April"}->[1] refers to (26, 27) 
#

# read in the files and store them in a data structure
#
for (my $i=0; $i<$n; $i++) {
    #print "Reading in set " . ($i+1) . " ...\n";
    open (LIST, "<$ARGV[$i]") or die ("Could not open $ARGV[$i]");
    my %set = ();
    while (<LIST>){
	chomp;
       	my @f = split;
	my $intval;
	if (!$simple && @f < 3) {
	    print "Simple mode with just lists of names.\n";
	    $simple = 1;
	}
	if ($simple) {
	    $intval = [1,1];
	} else {
	    $intval = [$f[1], $f[2]];
	}
	if (exists $set{$f[0]}){
	    push @{$set{$f[0]}}, $intval;
	} else {
	    $set{$f[0]} = [$intval];
	}
	$alllabels{$f[0]} = 0 if (!defined($alllabels{$f[0]}));
	$alllabels{$f[0]}++;
    }
    $rawsets[$i] = \%set;
    close LIST;
    $sets[$i] = {};
}


# process the raw data and create a new sets
# where complete labeled segments are the identifiers
#
# set1     *************
# set2           ***********        ******
# lbounds  |     |      |   |       |     |
#
#
foreach my $label (keys %alllabels) {
    # print "Finding boundaries for $label ";
    my %lbounds = ();
    # construct lbounds by combining all boundaries from all sets
    for (my $i=0; $i<$n; $i++) {
	foreach my $item (@{$rawsets[$i]->{$label}}) {
	    $lbounds{$item->[0]}++;
	    $lbounds{$item->[1]+1}++;
	}
    }
    my @slbounds = sort { $a <=> $b } keys %lbounds;
    # print join ", ", @slbounds;
    # print "\n";
    # now split up intervals at boundaries
    for (my $i=0; $i<$n; $i++) {
	# print "set$i ";
	foreach my $item (@{$rawsets[$i]->{$label}}) {
	    my $start = $item->[0];
	    my $end = $item->[1];
	    my $j = findIdx(\@slbounds, $start);
	    while($slbounds[$j] <= $end){
		my $key = $label . "\t" . $slbounds[$j] . "\t" . ($slbounds[$j+1]-1);
		$sets[$i]->{$key}++;
		$union{$key}++;
		$j++;
	    }
	}
	# print "\n";
    }
}

## individual list sizes
print "individual set sizes:\n";
for (my $i=0; $i<$n; $i++) {
    printf "set " . ($i+1) . " ";
    printf ("%30s : ", $ARGV[$i]);
    my $size = 0;
    foreach my $item (keys %{$sets[$i]}) {
	$size += weight($item);
    }
    print "$size\n";
}
printf "      ";
printf ("%30s : ", "total (union of sets)");
my $size = 0;
foreach my $item (keys %union) {
    $size += weight($item);
}
print "$size\n";

## power combinations
print "\npower combinations\n";
print "set         ";
for (my $i=0; $i<$n; $i++){
    printf ("%2d", $i+1);
}
print " | # of elements\n----------------------------------\n";
for (my $j=1; $j < (2**$n); $j++) {
   
    my @elem;
    my $combstr;
    print "combination ";
    for (my $i=0; $i<$n; $i++){
	$elem[$i] = ($j & 2**$i)? 1:0;
	$combstr .= sprintf ("%2d", $elem[$i]); 
    }
    print "$combstr | ";
    $combstr =~ s/ //g;
    if ($outfiles) {
	open (OUT, ">combset.$combstr.lst") or die ("Cound not open combset.$combstr.lst for writing.");
    }
    # now count the number of elements in the combination given by @elem
    my $k=0;
    foreach my $item (keys %union) {
	my $is_element = 1;
	for (my $i=0; $i<$n && $is_element; $i++){
	    my $in_set = ($sets[$i]->{$item}>0);
	    if (! ($in_set == $elem[$i])){
		$is_element = 0;
	    }
	}
	if ($is_element) {
	    if ($outfiles) {
		if ($simple) {
		    my @f = split "\t", $item;
		    print OUT "$f[0]\n";
		} else {
		    print OUT "$item\n";
		}
	    }
	    $k += weight($item);
	}
    }
    if ($outfiles) {
	close OUT;
    }
    print $k . "\n";
}

print "\n\nExplanation:\n";
print "1: is in the set\n";
print "0: is not in the set\n";
print "For example, if you have n=3 sets, then the line\n";
print "combination 1 0 0\n";
print "Shows the number of items that are in set 1 but neither in set 2 nor in set 3.\n";
print "All 2^n-1 = ". (2**$n-1). " combinations add up to the total number of elements.\n";



sub findIdx {
    my $bounds = shift;
    my $elem = shift;
    # binary search for efficiency
    my $a=0;
    my $b=@{$bounds}-1;
    while ($b>$a) {
	return $a if ($bounds->[$a] == $elem);
	my $c = int(($a+$b+1)/2);
	if ($bounds->[$c] <= $elem) {
	    $a=$c;
	} else {
	    $b=$c;
	}
    }
    print STDERR "Error: Index not found\n";
}

sub weight {
    my $key = shift;
    if ($key =~ /\t(\d+)\t(\d+)$/) {
	return $2-$1+1;
    } else {
	print STDERR "Error in weight function: key=$key\n";
    }
}
