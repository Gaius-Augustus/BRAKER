#!/usr/bin/perl

##############################################################
# filterGenes
# filter genes from a genbank flat file database
# usage: fileterGenes namefile dbfile
#
#
# Mario Stanke & Katharina Hoff, last modification Feb 19 2018
##############################################################

if ($#ARGV != 1) {
    print "usage: filterGenesIn.pl namefile dbfile\n";
    print "names of the loci to be kept come from\n";
    print "the first parameter. Only the the first of identical loci is kept\n";
    exit;
}
$origfilename = $ARGV[1];
$goodfilename = $ARGV[0];
open(goodfile, "<$goodfilename") || die "ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCouldn't open name file";

while(<goodfile>){
   /.*/;
   $goodlist{$&}=1;
}


open(origfile, "<$origfilename") || die "ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCouldn't open dbfile\n";


$/="\n//\n";
while(<origfile>) {
    $gendaten=$_;
    m/^LOCUS +(\S+) .*/;
    $genname=$1;
    if (exists($goodlist{$genname})) {
	print "$gendaten";
	delete($goodlist{$genname});
    }
}
