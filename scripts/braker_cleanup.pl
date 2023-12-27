#!/usr/bin/env perl

####################################################################################################
#                                                                                                  #
# braker_cleanup.pl                                                                                #
# Script that deletes all files produced by BRAKER that are not deleted by the braker.pl script,   #
# itself, and that are usually not required for downstream analysis, unless you wish to debug      #
# something.                                                                                       #
#                                                                                                  #
# Author: Katharina Hoff                                                                           #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
####################################################################################################

use Getopt::Long;
use File::Path;
use File::Spec::Functions qw(rel2abs);
use strict;
use warnings;
use File::Copy;

my $usage = <<'ENDUSAGE';

braker_cleanup.pl     delete braker.pl output files that are usually not required for 
                      downstream analysis

SYNOPSIS

braker_cleanup.pl --wdir=WDIR

OPTIONS

--wdir=WDIR         output directory of braker.pl job

--help	            Display this help message

ENDUSAGE

my $help;
my $wdir;
my $homedir = $ENV{HOME};

GetOptions(
	'wdir=s'  => \$wdir,
    'help!' => \$help
);

if ($help) {
	print $usage;
    exit(0);
}

if(not(defined($wdir))){
	print "ERROR: in file " . __FILE__ ." at line "
            . __LINE__ . "\n" . "No Working directory provided!"
            . "(option --wdir=WDIR)\n";
    print $usage;
	exit(1)
}

if($wdir =~ m/^~/){
	$wdir =~ s/^~//;
	$wdir = File::Spec->catfile($homedir, $wdir);
}

if(not($wdir =~ m/^\//)){
	$wdir = rel2abs($wdir);
}

if(not(-d $wdir)){
	my @t = split(/\//, $wdir);
	$wdir = "";
	for(my $i=0; $i<(scalar(@t) - 1); $i++){
		$wdir .= "/".$t[$i]
	}
	if(not(-d $wdir)){
		print "ERROR: in file " . __FILE__ ." at line "
            . __LINE__ . "\n" . "wdir $wdir is not a directory!\n";
	exit(1)
	}
}

my @files = ("firsttest.stdout", "genome.fa", "getAnnoFasta.augustus.ab_initio.stdout", 
	"getAnnoFasta.augustus.hints.stdout", "optimize_augustus.stdout", "secondtest.stdout",
	"train.gb", "GeneMark-ES/genemark.d.gtf", "GeneMark-ES/genemark.f.good.gtf", 
	"GeneMark-ES/gmes.log", "GeneMark-ES/logfile", "GeneMark-ES/run.cfg", 
	"GeneMark-ES/training.fna", "GeneMark-ES.stdout", "GeneMark-ET/genemark.d.gtf", 
	"GeneMark-ET/genemark.f.good.gtf", "GeneMark-ET/genemark.c.gtf", "GeneMark-ET/genemark.average_gene_length.out",
	"GeneMark-ET/gmes.log", "GeneMark-ET/logfile", "GeneMark-ET/run.cfg", 
	"GeneMark-ET/training.fna", "GeneMark-ET/genemark.f.bad.gtf", "GeneMark-ET.stdout", "filterGenemark.stdout", 
	"GeneMark-EP/genemark.d.gtf", 
	"GeneMark-EP/genemark.f.good.gtf", "GeneMark-EP/genemark.c.gtf", "GeneMark-EP/genemark.average_gene_length.out",
	"GeneMark-EP/gmes.log", "GeneMark-EP/logfile", "GeneMark-EP/run.cfg", 
	"GeneMark-EP/training.fna", "GeneMark-EP/genemark.f.bad.gtf", "GeneMark-EP.stdout",
	"GeneMark-ETP/genemark.d.gtf", 
	"GeneMark-ETP/genemark.f.good.gtf", "GeneMark-ETP/genemark.c.gtf", "GeneMark-ETP/genemark.average_gene_length.out",
	"GeneMark-ETP/gmes.log", "GeneMark-ETP/logfile", "GeneMark-ETP/run.cfg", "GeneMark-ETP/proteins.fa/model/set.out",
    "GeneMark-ETP/proteins.fa/model/parse.es_c", "GeneMark-ETP/proteins.fa/model/non.seq", "GeneMark-ETP/proteins.fa/nonhc/nonhc.fasta",
    "GeneMark-ETP/rnaseq/hints/proteins.fa/prothint/nuc.fasta", 
	"GeneMark-ETP/training.fna", "GeneMark-ETP/genemark.f.bad.gtf", "GeneMark-ETP.stdout",
	"aug_hints.lst", "aa2nonred.stdout", "augustus.hints.tmp.gtf", "bothutr.lst", "fourthtest.stdout", "gbFilterEtraining.stdout",
	"genes.gtf", "genes_in_gb.gtf", "gff2gbSmallDNA.utr.stdout", "hints.job.lst", "merged.bam", "merged.s.bam", "merged.wig", 
	"randomSplit_utr1.log", "randomSplit_utr2.log", "rnaseq2utr.stdout", "rnaseq.utr.hints", "samtools_sort_before_wig.stdout",
	"singlecds.hints", "stops.and.starts.gff", "train.gb.test", "train.gb.train", "train.gb.train.test",
	"train.gb.train.train", "traingenes.good.fa", "utr.aa2nonred.stdout",
	"utr.gb", "utr.gb.test", "utr.gb.train", "utr.gb.train.test", "utr.gb.train.train", "utr_genes_in_gb.fa",
	"utr_genes_in_gb.nr.fa", "utr.nonred.loci.lst", "utrs.gff", "downsample_traingenes.log", "firstetraining.stdout", 
	"secondetraining.stdout", "startAlign_gth.log", "protein_alignment_gth.gff3", "ex1.cfg", "getAnnoFastaFromJoingenes.augustus.hints.stdout", 
	"genome.fa.cidx", "getAnnoFastaFromJoingenes.augustus.hints_tmp.stdout", "getAnnoFastaFromJoingenes.augustus.ab_initio_tmp.stdout",
	"augustus.ab_initio.tmp.gtf", "augustus.ab_initio.gff", "augustus.hints.gff", "augustus.hints.tmp.gtf",
	"getAnnoFastaFromJoingenes.augustus.hints_hints.stdout", "getAnnoFastaFromJoingenes.augustus.ab_initio_.stdout",
	"getAnnoFastaFromJoingenes.augustus.hints_.stdout", "startAlign.stdout", "augustus.hints_iter1.aa", "augustus.hints_iter1.codingseq",
	"augustus.hints_iter1.gff", "augustus.hints_iter1.gtf", "augustus.tmp1.gff", "augustus.tmp2.gff", "cmd.log", "etrain.bad.lst", "gene_stat.yaml",
	"good_genes.lst", "hintsfile_iter1.gff", "nonred.loci.lst", "nuc.fasta", "prevHints.gff", "proteins.fa", "seed_proteins.faa",
    "compleasm_hints.gff", "train.f.gb", "traingenes.good.nr.fa", "uniqueSeeds.gtf", "braker.gtf_temp", "genemark_hintsfile.gff", "gc_content.out", "evidence.gff", "genemark_evidence.gff");

my @dirs = ("GeneMark-ES/data", "GeneMark-ES/info", "GeneMark-ES/output", "GeneMark-ES/run",
	"GeneMark-ET/data", "GeneMark-ET/info", "GeneMark-ET/output", "GeneMark-ET/run",
	"GeneMark-EP/data", "GeneMark-EP/info", "GeneMark-EP/output", "GeneMark-EP/run",
	"GeneMark-ETP/etp_data", "GeneMark-ETP/rnaseq/gmst/", "GeneMark-ETP/rnaseq/hisat2/", 
 #   "GeneMark-ETP/rnaseq/stringtie/", 
	"GeneMark-ETP/arx", "GeneMark-ETP/proteins.fa/hc/regions/",
	"GeneMark-ETP/proteins.fa/nonhc/prothint","GeneMark-ETP/rnaseq/hints/proteins.fa/prothint/Spaln",
    "GeneMark-ETP/proteins.fa/nonhc/for_prothint", "GeneMark-ETP/proteins.fa/nonhc/pred_m",
    "GeneMark-ETP/proteins.fa/nonhc/regions", "GeneMark-ETP/proteins.fa/etr", "GeneMark-ETP/proteins.fa/nonhc/hc", 
    "GeneMark-ETP/rnaseq/hints/proteins.fa/prothint/diamond", "rnaseq/", "GeneMark-ETP/data",
    "genome_split", "Spaln", "augustus_tmp", "augustus_files_hints", "diamond", "mb_downloads");

foreach(@files){
	if(-e $wdir."/".$_){
		print "Deleting file ".$wdir."/".$_."\n";
		unlink($wdir."/".$_);
	}
}

foreach(@dirs){
	if(-d $wdir."/".$_){
		print "Deleting directory ".$wdir."/".$_."\n";
		rmtree($wdir."/".$_);
	}
}

my $gth_index = 0;
while(-d $wdir."/align_gth".$gth_index){
	print "Deleting directory ".$wdir."/align_gth".$gth_index."\n";
	rmtree($wdir."/align_gth".$gth_index);
	$gth_index = $gth_index + 1;
}

# move augustus.* output files into a subdirectory
opendir my $dir, "$wdir" or die "Cannot open directory: $wdir!";
@files = readdir $dir;
closedir $dir;

mkdir("$wdir/Augustus", 0700) or die("Failed to create directory Augustus!\n");

foreach(@files){
	if($_ =~ m/augustus\./){
		move($wdir."/".$_, $wdir."/Augustus/".$_) or die("Failed to move file $wdir/$_ to $wdir/Augustus/$_!\n");
	}
}
