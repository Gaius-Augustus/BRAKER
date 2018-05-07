#!/usr/bin/env perl

use Getopt::Long;
use File::Compare;
use File::Path qw(make_path rmtree);
use Module::Load::Conditional qw(can_load check_install requires);
use Scalar::Util::Numeric qw(isint);
use POSIX qw(floor);
use List::Util qw[min max];
use Parallel::ForkManager;
use FindBin;
use lib "$FindBin::RealBin/.";
use File::Which;                    # exports which()
use File::Which qw(which where);    # exports which() and where()
use Cwd;
use Cwd 'abs_path';
use File::Spec::Functions qw(rel2abs);
use File::Basename qw(dirname basename);
use File::Copy;
use helpMod
    qw( find checkFile formatDetector relToAbs setParInConfig uptodate gtf2fasta clean_abort );
use Term::ANSIColor qw(:constants);
use strict;
use warnings;
############## CODE TO SIMULATE BRAKER RUN #####################################

my $otherfilesDir = "/home/katharina/utr_test";
my $errorfilesDir = $otherfilesDir;
my $AUGUSTUS_CONFIG_PATH = "/home/katharina/SVN/augustus/trunks/config";
my $hintsfile = $otherfilesDir."/hintsfile.gff";
my $genome = $otherfilesDir."/genome.fa";
my $BAMTOOLS_BIN_PATH = "/home/katharina/git/bamtools/build/src/toolkit/";
my $rnaseq2utrPath = "/home/katharina/SVN/utrrnaseq/trunks/Debug/utrrnaseq";# full path with toolname
my $bam2wigPath = "/home/katharina/SVN/augustus/trunks/auxprogs/bam2wig/bam2wig";#
my $v = 4;
my $species = "Sp_1";
my @bam = ("/home/katharina/utr_test/RNAseq.bam");
my $cmdString;
my $perlCmdString;
my $nice = 1;
my $rnaseq2utr_args = "";
my $string;
my $AUGUSTUS_BIN_PATH = "/home/katharina/SVN/augustus/trunks/bin/";
my $AUGUSTUS_SCRIPTS_PATH = "/home/katharina/SVN/augustus/trunks/scripts/";
my $flanking_DNA = 1000; # TODO: make sure that this is global in braker.pl
my $rounds = 3;
my $CPU = 8;
my $BLAST_PATH = "/usr/bin/"
open(LOG, ">", "/home/katharina/utr.log");
train_utr();



####################### train_utr ##############################################
# * train UTR parameters for AUGUSTUS
################################################################################

sub train_utr {
    print LOG "\# " . (localtime) . "Training AUGUSTUS UTR parameters\n"
        if ( $v > 2 );
    print LOG "\# "
        . (localtime)
        . ": Move augustus predictions to *.noUTR.* files prior UTR training:\n"
        if ( $v > 3 );
    my $loci = 0;
    my $testSetSize = 0;
    my $onlyTrainSize = 0;
    # store predictions without UTRs, revert later if UTR model does not improve
    # predictions
    print LOG "mv $otherfilesDir/augustus.hints.gff "
        . "$otherfilesDir/augustus.hints.noUtr.gff\n" if ( $v > 3 );
    move( "$otherfilesDir/augustus.hints.gff",
        "$otherfilesDir/augustus.hints.noUtr.gff" );
    print LOG "mv $otherfilesDir/augustus.hints.gtf "
        . "$otherfilesDir/augustus.hints.noUtr.gtf\n" if ( $v > 3 );
    move( "$otherfilesDir/augustus.hints.gtf",
        "$otherfilesDir/augustus.hints.noUtr.gtf" );
    if(-e "$otherfilesDir/augustus.hints.aa"){
        print LOG "mv $otherfilesDir/augustus.hints.aa "
            . "$otherfilesDir/augustus.hints.noUtr.aa\n" if ( $v > 3 );
        move( "$otherfilesDir/augustus.hints.aa",
            "$otherfilesDir/augustus.hints.noUtr.aa" );
    }
    if(-e "$otherfilesDir/augustus.hints.codingseq"){
        print LOG "mv $otherfilesDir/augustus.hints.codingseq "
            . "$otherfilesDir/augustus.hints.noUtr.codingseq\n" if ( $v > 3 );
        move( "$otherfilesDir/augustus.hints.codingseq",
            "$otherfilesDir/augustus.hints.noUtr.codingseq" );
    }
    if(-e "$otherfilesDir/augustus.hints.cdsexons"){
        print LOG "mv $otherfilesDir/augustus.hints.cdsexons "
            . "$otherfilesDir/augustus.hints.noUtr.cdsexons\n" if ( $v > 3 );
        move( "$otherfilesDir/augustus.hints.cdsexons",
            "$otherfilesDir/augustus.hints.noUtr.cdsexons" );
    }

    # copy species parameter files, revert later if UTR model does not improve
    # predictions
    print LOG "\# " . (localtime)
        . ": Create backup of current species parameters:\n" if ( $v > 3 );
    for ( ( "$species" . "_exon_probs.pbl",
            "$species" . "_igenic_probs.pbl",
            "$species" . "_intron_probs.pbl" ) ) {
        print LOG "cp $AUGUSTUS_CONFIG_PATH/species/$species/$_ "
            . "$AUGUSTUS_CONFIG_PATH/species/$species/$_.noUTR\n" if ( $v > 3 );
        copy( "$AUGUSTUS_CONFIG_PATH/species/$species/$_",
              "$AUGUSTUS_CONFIG_PATH/species/$species/$_.noUTR" )
            or die( "ERROR in file " . __FILE__ . " at line " . __LINE__
                . "\nCopy of $AUGUSTUS_CONFIG_PATH/species/$species/$_ to "
                . "$AUGUSTUS_CONFIG_PATH/species/$species/$_.noUTR failed!\n" );
    }
    chdir($otherfilesDir) or die( "ERROR in file " . __FILE__ . " at line "
            . __LINE__
            . "\nCould not change into directory $otherfilesDir!\n" );

    # search all start and stop codons from augustus.noUtr.gtf and write them
    # to the file stops.and.starts.gff
    if ( !uptodate( ["$otherfilesDir/augustus.hints.noUtr.gtf"],
        ["$otherfilesDir/stops.and.starts.gff"] ) ) {
        print LOG "\# " . (localtime)
            . ": extracting all stop and start codons from "
            . "augustus.hints.noUtr.gtf to stops.and.starts.gff\n"
            if ( $v > 3 );
        my %nonRedundantCodons;
        my @tmpGffLine;
        open( AUG, "<", "$otherfilesDir/augustus.hints.noUtr.gtf" )
            or die( "ERROR in file " . __FILE__ . " at line " . __LINE__
                . "\nCould not open file "
                . "$otherfilesDir/augustus.hints.noUtr.gtf!\n" );
        while ( defined( my $i = <AUG> ) ) {
       # TODO: we are not dealing with redundancy, correctly. Discarding
       #       duplicates is not the optimal solution because later, we filter
       #       for genes where both codons have UTR models. However, at this
       #       point in time, we ignore this matter and hope that we are left
       #       with a sufficient number of training examples.
            if ( $i =~ /\t(start_codon|stop_codon)\t/ ) {
                @tmpGffLine = split( /\t/, $i );
                if ( not( defined( $nonRedundantCodons{"$tmpGffLine[0]_$tmpGffLine[3]_$tmpGffLine[4]_$tmpGffLine[6]"} ) ) ) {
                    $nonRedundantCodons{"$tmpGffLine[0]_$tmpGffLine[3]_$tmpGffLine[4]_$tmpGffLine[6]" } = $i;
                }
            }
        }
        close(AUG) or die( "ERROR in file " . __FILE__ . " at line " . __LINE__
            . "\nCould not close file "
            . "$otherfilesDir/augustus.hints.noUtr.gtf!\n" );
        open( CODON, ">", "$otherfilesDir/stops.and.starts.gff" )
            or die( "ERROR in file " . __FILE__ . " at line " . __LINE__
                . "\nCould not open file "
                . "$otherfilesDir/stops.and.starts.gff!\n" );
        foreach my $key ( keys %nonRedundantCodons ) {
            print CODON $nonRedundantCodons{$key};
        }
        close(CODON) or die( "ERROR in file " . __FILE__ . " at line "
            . __LINE__ . "\nCould not close file "
            . "$otherfilesDir/stops.and.starts.gff!\n" );
    }

    if ( !uptodate( ["$hintsfile"], ["$otherfilesDir/rnaseq.utr.hints"] ) ) {
      # TODO: currently, only using AT-AG, not AC-AG or any other splice site.
      #       Possibly extend to other splice patterns.
        print LOG "\# " . (localtime)
            . ": filtering RNA-Seq hints for valid splice site AT-AG, storing "
            . "in $otherfilesDir/rnsaeq.utr.hints\n" if ( $v > 3 );
        #TODO: do not load entire genome in memory, process chromsome-wise!
        my %genome_hash;
        my $hash_key;
        open( FASTA, "<", $genome ) or die( "ERROR in file " . __FILE__
            . " at line " . __LINE__ . "\nCould not open file $genome!\n" );
        LINE: while ( my $line = <FASTA> ) {
            next LINE if $line =~ m/^#/;
            if ( $line =~ /^>/ ) {
                chomp($line);
                $hash_key = substr( $line, 1, length($line) - 1 );
            } else {
                $line =~ s/[\x0A\x0D]+//g;
                $line =~ s/(\s+)(\n)(\r)//g;
                $line = uc($line);
                $genome_hash{$hash_key} .= $line;
            }
        }
        close(FASTA) or die( "ERROR in file " . __FILE__ . " at line "
                . __LINE__ . "\nCould not close file $genome!\n" );
        open( HINTS, "<", $hintsfile ) or die( "ERROR in file " . __FILE__
            . " at line " . __LINE__ . "\nCould not open file $hintsfile!\n" );
        my @gff;
        my ( $siteA, $siteB, $given, $lastCol );
        my $splice = "GTAG";
        open( UTRHINTS, ">", "$otherfilesDir/rnaseq.utr.hints" )
            or die( "ERROR in file " . __FILE__ . " at line " . __LINE__
            . "\nCould not open file rnaseq.utr.hints!\n" );
        while ( my $line = <HINTS> ) {
            @gff = split( /\t/, $line );
            if ( ( $gff[1] eq "b2h" ) && ( $gff[2] eq "intron" ) ){
                $siteA = substr( $genome_hash{ $gff[0] }, ( $gff[3] - 1 ), 2 );
                $siteB = substr( $genome_hash{ $gff[0] }, ( $gff[4] - 2 ), 2 );
                $given = $siteA . $siteB;
                if ( $gff[8] =~ m/mult=(\d+)/ ) {
                    $lastCol = "mult=".$1."_$splice\n";
                }else {
                    $lastCol = "mult=1_$splice\n";
                }
                if ( uc($given) =~ m/$splice/ ) {
                    for(my $i = 0; $i < 8; $i++){
                        if($i != 6 ) {
                            print UTRHINTS $gff[$i] . "\t";
                        }else{
                            print UTRHINTS "+\t";
                        }
                    }
                    print UTRHINTS $lastCol;
                } else {
                    $given = reverse $given;
                    $given =~ tr/ACGTacgt/TGCAtgca/;
                    if ( uc($given) =~ m/$splice/ ) {
                        for(my $i = 0; $i < 8; $i++){
                        if($i != 6 ) {
                            print UTRHINTS $gff[$i] . "\t";
                        }else{
                            print UTRHINTS "-\t";
                        }
                    }
                    print UTRHINTS $lastCol;
                    }
                }
            }
        }
        close(UTRHINTS) or die( "ERROR in file " . __FILE__ . " at line "
            . __LINE__ . "\nCould not close file "
            . "$otherfilesDir/rnaseq.utr.hints!\n" );
        close(HINTS) or die( "ERROR in file " . __FILE__ . " at line "
            . __LINE__ . "\nCould not close file $hintsfile!\n" );
    }

    # create wiggle file from bam files
    if ( !uptodate( ["$hintsfile"], ["$otherfilesDir/merged.wig"] ) ) {
        print LOG "\# " . (localtime)
                . ": Converting bam files to wiggle file merged.wig\n"
                if ( $v > 3 );
        if ( scalar(@bam) > 1 ) {
            print LOG "\# " . (localtime)
                . ": For conversion, merge multiple bam files\n" if ( $v > 3 );
            $cmdString = "";
            if ($nice) {
                $cmdString .= "nice ";
            }
            $cmdString .= "$BAMTOOLS_BIN_PATH/bamtools merge ";
            foreach (@bam) {
                chomp;
                $cmdString .= "-in $_ ";
            }
            $cmdString .= "-out $otherfilesDir/merged.bam "
                       .  "1> $otherfilesDir/bam.merge.log "
                       .  "2> $errorfilesDir/bam.merge.err";
            print LOG "\n$cmdString\n\n" if ( $v > 3 );
            system("$cmdString") or die( "ERROR in file " . __FILE__
                    . " at line " . __LINE__
                    . "\nFailed to execute: $cmdString!\n" );
        } else {
            print LOG "\# " . (localtime) . ":  For conversion, creating "
                . "softlink to bam file $bam[0]...\n" if ( $v > 3 );
            $cmdString = "ln -s $bam[0] $otherfilesDir/merged.bam";
            print LOG "$cmdString\n" if ( $v > 3 );
            system($cmdString) == 0 or die( "ERROR in file " . __FILE__
                . " at line " . __LINE__
                . "\nFailed to execute: $cmdString!\n" );
        }
        print LOG "\# " . (localtime) . ": Creating wiggle file...\n"
            if ( $v > 3 );
        $cmdString = "";
        if ($nice) {
            $cmdString .= "nice ";
        }
        $cmdString .= "$bam2wigPath merged.bam >$otherfilesDir/merged.wig "
                   .   "2> $errorfilesDir/bam2wig.err";
        print LOG "\n$cmdString\n" if ( $v > 3 );
        system("$cmdString") == 0 or die( "ERROR in file " . __FILE__
            . " at line " . __LINE__ . "\nFailed to execute: $cmdString!\n" );
    }

    # call utrrnaseq
    if ( !uptodate( ["$hintsfile"], ["$otherfilesDir/utrs.gff"] ) ) {
        print LOG "\# " . (localtime) . ": Creating $otherfilesDir/utrs.gff\n"
            if ( $v > 3 );
        $cmdString = "";
        if ($nice) {
            $cmdString .= "nice ";
        }
        $cmdString .= "$rnaseq2utrPath --in-scaffold-file $genome "
                   .  "-C $otherfilesDir/stops.and.starts.gff "
                   .  "-I $otherfilesDir/rnaseq.utr.hints "
                   .  "-W $otherfilesDir/merged.wig "
                   .  "-o $otherfilesDir/utrs.gff "
                   .  "$rnaseq2utr_args 1> $otherfilesDir/rnaseq2utr.log "
                   .  "2> $errorfilesDir/rnaseq2utr.err";
        print LOG "\n$cmdString\n" if ( $v > 3 );
        system("$cmdString") == 0 or die( "ERROR in file " . __FILE__
            . " at line " . __LINE__ . "\nFailed to execute: $cmdString!\n" );
    }

    # create genbank file with genes that have two utrs
    if (!uptodate( [ "$otherfilesDir/utrs.gff",    "$otherfilesDir/augustus.hints.noUtr.gtf" ],
            [ "$otherfilesDir/bothutr.lst", "$otherfilesDir/bothutr.test.gb" ] ) ) {
        print LOG "\# " . (localtime) . ": Creating gb file for UTR training\n"
            if ( $v > 3 );

        # extract subset of genes, where we have both UTRs
        open( UTR, "<", "$otherfilesDir/utrs.gff" ) or die( "ERROR in file "
            . __FILE__ . " at line " . __LINE__
            . "\nCan not open file $otherfilesDir/utrs.gff!\n" );
        my %both;
        while (<UTR>) {
            my @t = split(/\t/);
            $t[8]=~m/transcript_id \"(\S+)\"/;
            my $txid = $1;
            if($t[2] =~ m/UTR/){
                $both{$txid}{$t[2]} = $_;
            }
        }
        close(UTR) or die( "ERROR in file " . __FILE__ . " at line " . __LINE__
            . "\nCould not close file $otherfilesDir/utrs.gff!\n" );

        open( BOTH, ">", "$otherfilesDir/bothutr.lst" ) or die( "ERROR in file "
            . __FILE__ . " at line " . __LINE__
            . "\nCan not open file $otherfilesDir/bothutr.lst!\n" );
        foreach(keys %both){
            if(defined($both{$_}{"3'-UTR"}) && defined($both{$_}{"5'-UTR"})){
                print BOTH $_."\n";
            }
        }
        close(BOTH) or die( "ERROR in file " . __FILE__ . " at line "
            . __LINE__
            . "\nCould not close file $otherfilesDir/bothutr.lst!\n" );
        $cmdString = "";
        if ($nice) {
            $cmdString .= "nice ";
        }
        $cmdString .= "cat $otherfilesDir/utrs.gff "
                   .  "$otherfilesDir/augustus.hints.noUtr.gtf | "
                   .  "grep -P \"(CDS|5'-UTR|3'-UTR)\" | "
                   .  "sort -n -k 4,4 | "
                   .  "sort -s -k 10,10 | sort -s -k 1,1 >"
                   .  "> $otherfilesDir/genes.gtf "
                   .  "2> $errorfilesDir/cat_utrs_augustus_noUtrs.err";
        print LOG "\n$cmdString\n" if ( $v > 3 );
        system($cmdString) == 0 or die( "ERROR in file " . __FILE__
            . " at line " . __LINE__
            . "\nFailed not execute $cmdString!\n" );

        $string = find(
            "gff2gbSmallDNA.pl",    $AUGUSTUS_BIN_PATH,
            $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH
        );
        $perlCmdString = "";
        if ($nice) {
            $perlCmdString .= "nice ";
        }
        $perlCmdString .= "perl $string $otherfilesDir/genes.gtf $genome "
                       .  "$flanking_DNA $otherfilesDir/utr.gb "
                       .  "--good=$otherfilesDir/bothutr.lst "
                       .  "1> $otherfilesDir/gff2gbSmallDNA.utr.stdout "
                       . "2> $errorfilesDir/gff2gbSmallDNA.utr.stderr";
        print LOG "\n$perlCmdString\n" if ( $v > 3 );
        system("$perlCmdString") == 0 or die( "ERROR in file " . __FILE__
            . " at line " . __LINE__
            . "\nFailed to execute: $perlCmdString!\n" );

        # assign LOCI locations to txIDs
        open (TRAINUTRGB1, "<", "$otherfilesDir/utr.gb") or die( 
            "ERROR in file " . __FILE__ . " at line " . __LINE__
            . "\nCould not open file $otherfilesDir/utr.gb!\n" );

        my %txInUtrGb1;
        my $txLocus;;
        while( <TRAINUTRGB1> ) {
            if ( $_ =~ m/LOCUS\s+(\S+)\s/ ) {
                $txLocus = $1;
            }elsif ( $_ =~ m/\/gene=\"(\S+)\"/ ) {
                $txInUtrGb1{$1} = $txLocus;
            }
        }
        close (TRAINUTRGB1) or die( 
            "ERROR in file " . __FILE__ . " at line " . __LINE__
            . "\nCould not close file $otherfilesDir/utr.gb!\n" );

        # find those genes in gtf that ended up in gb file
        open(GENES, "<", "$otherfilesDir/genes.gtf") or die( "ERROR in file " 
            . __FILE__ . " at line " . __LINE__
            . "\nCould not open file $otherfilesDir/genes.lst!\n" );
        open(GENESINGB, ">", "$otherfilesDir/genes_in_gb.gtf") or die( 
            "ERROR in file " . __FILE__ . " at line " . __LINE__
            . "\nCould not open file $otherfilesDir/genes_in_gb.gtf!\n" );
        while(<GENES>){
            $_=~m/transcript_id \"(\S+")\"/;
            if(defined($txInUtrGb1{$1})){
                print GENESINGB $_;
            }
        }
        close(GENESINGB) or die( "ERROR in file " 
            . __FILE__ . " at line " . __LINE__
            . "\nCould not close file $otherfilesDir/genes_in_gb.gtf!\n" );
        close(GENES) or die( "ERROR in file " 
            . __FILE__ . " at line " . __LINE__
            . "\nCould not close file $otherfilesDir/genes.lst!\n" );

        # convert those training genes to protein fasta file
        gtf2fasta ($genome, "$otherfilesDir/genes_in_gb.gtf",
            "$otherfilesDir/utr_genes_in_gb.fa");

        # blast good training genes to exclude redundant sequences
        $string = find(
            "aa2nonred.pl",       $AUGUSTUS_BIN_PATH,
            $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH
        );
        $perlCmdString = "";
        if ($nice) {
            $perlCmdString .= "nice ";
        }
        $perlCmdString .= "perl $string $otherfilesDir/utr_genes_in_gb.fa "
                       .  "$otherfilesDir/utr_genes_in_gb.nr.fa "
                       .  "--BLAST_PATH=$BLAST_PATH --cores=$CPU "
                       .  "1> $otherfilesDir/utr.aa2nonred.out "
                       .  "2> $errorfilesDir/utr.aa2nonred.stderr";
        print LOG "\# "
            . (localtime)
            . ": BLAST training gene structures (with UTRs) against "
            . "themselves:\n" if ($v > 3);
        print LOG "$perlCmdString\n\n" if ($v > 3);
        system("$perlCmdString") == 0
            or die( "ERROR in file " . __FILE__
            . " at line " . __LINE__
            . "\nFailed to execute: $perlCmdString!\n" );

        # parse output of blast
        my %nonRed;
        open (BLASTOUT, "<", "$otherfilesDir/utr_genes_in_gb.nr.fa") or
             die( "ERROR in file "
            . __FILE__ . " at line " . __LINE__
            . "\nCould not open file $otherfilesDir/utr_genes_in_gb.nr.fa!\n" );
        while ( <BLASTOUT> ) {
            chomp;
            if($_ =~ m/^\>(\S+)/){
                $nonRed{$1} = 1;
            }
        }
        close (BLASTOUT) or  die( "ERROR in file "
            . __FILE__ . " at line " . __LINE__
            . "\nCould not close file $otherfilesDir/utr_genes_in_gb.nr.fa!\n" );

        open ( NONREDLOCI, ">", "$otherfilesDir/utr.nonred.loci.lst") or
            die( "ERROR in file "
            . __FILE__ . " at line " . __LINE__
            . "\nCould not open file $otherfilesDir/utr.nonred.loci.lst!\n" );
        foreach ( keys %nonRed ) {
            print NONREDLOCI $txInUtrGb1{$_}."\n";
        }
        close (NONREDLOCI) or  die( "ERROR in file "
            . __FILE__ . " at line " . __LINE__
            . "\nCould not close file $otherfilesDir/utr.nonred.loci.lst!\n" );

        # filter utr.gb file for nonredundant loci
        $string = find(
            "filterGenesIn.pl",       $AUGUSTUS_BIN_PATH,
            $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH
        );
        $errorfile     = "$errorfilesDir/utr.filterGenesIn.stderr";
        $perlCmdString = "";
        if ($nice) {
            $perlCmdString .= "nice ";
        }
        $perlCmdString .= "perl $string $otherfilesDir/utr.nonred.loci.lst "
                       .  "$otherfilesDir/utr.gb 1> $otherfilesDir/utr.nr.gb "
                       .  "2>$errorfile";
        print LOG "\# " . (localtime)
            . ": Filtering nonredundant loci into $otherfilesDir/utr.nr.gb:\n" 
            if ($v > 3);
        print LOG "$perlCmdString\n\n" if ($v > 3);
        system("$perlCmdString") == 0
            or die( "ERROR in file " . __FILE__
            . " at line " . __LINE__
            . "\nFailed not execute $perlCmdString!\n" );

        # count how many genes are still in utr.nr.gb
        if($v > 3) {
            open (TRAINUTRGB2, "<", "$otherfilesDir/utr.nr.gb") or
                die( "ERROR in file " . __FILE__ . " at line " . __LINE__
                . "\nCould not open file $otherfilesDir/utr.nr.gb!\n" );
            my $nLociUtrGb2 = 0;
            while ( <TRAINUTRGB2> ) {
                if($_ =~ m/LOCUS/) {
                    $nLociUtrGb2++;
                }
            }
            close (TRAINUTRGB2) or  die( "ERROR in file " . __FILE__ 
                . " at line " . __LINE__
                . "\nCould not close file $otherfilesDir/utr.nr.gb!\n" );
            print LOG "\# "
                    . (localtime)
                    . ": $otherfilesDir/utr.nr.gb file contains $nLociUtrGb2 genes.\n";
        }

        # move nonredundant file utr.nr.gb to utr.gb
        $cmdString = "mv $otherfilesDir/utr.nr.gb $otherfilesDir/utr.gb";
        print LOG  "\# " . (localtime) . ": Moving utr.nr.gb to utr.gb file:\n"
            if($v>3);
        print LOG $cmdString."\n";
        system("$cmdString") == 0 or die( "ERROR in file " . __FILE__
            . " at line " . __LINE__
            . "\nFailed to execute: $cmdString!\n" );

        # create an utr.onlytrain.gb if more than 200 UTR training genes
        open(UTRGB, "<", "$otherfilesDir/utr.gb") or die( "ERROR in file "
            . __FILE__ . " at line " . __LINE__
            . "\nCould not open file $otherfilesDir/utr.gb!\n" );
        while(<UTRGB>){
            if($_ =~ m/LOCUS/){
                $loci++
            }
        }
        close(UTRGB) or die( "ERROR in file "
            . __FILE__ . " at line " . __LINE__
            . "\nCould not close file $otherfilesDir/utr.gb!\n" );
        if($loci < 50){
            die("ERROR in file " . __FILE__ . " at line " . __LINE__
            . "Number of UTR training loci is smaller than 50, aborting "
            . "UTR training! If this is the only error message, the "
            . "AUGUSTUS parameters for your species were optimized ok, "
            . "but you are lacking UTR parameters. Do not attempt to "
            . "predict genes with UTRs for this species using the current "
            . "parameter set!\n");
        }elsif( ($loci >= 50) && ($loci <= 1000) ){
            $testSetSize = floor($loci * 0.2);
            if(($loci - $testSetSize) > 200){
                $onlyTrainSize = 200;
            }
        }elsif($loci > 1000){
            $testSetSize = 200;
            $onlyTrainSize = 200;
        }

        # create hold out test set
        $string = find(
            "randomSplit.pl",       $AUGUSTUS_BIN_PATH,
            $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH
        );
        print LOG "Found script $string.\n" if ( $v > 3 );
        $perlCmdString = "perl $string $otherfilesDir/utr.gb $testSetSize "
                       . "1> $otherfilesDir/randomSplit_utr1.log "
                       . "2> $errorfilesDir/randomSplit_utr1.err";
        print LOG "\n$perlCmdString\n" if ( $v > 3 );
        system("$perlCmdString") == 0 or die( "ERROR in file " . __FILE__
            . " at line " . __LINE__
            . "\nFailed to execute: $perlCmdString!\n" );

        # create onlytrain if applicable
        if($onlyTrainSize > 0){
             $string = find(
                "randomSplit.pl",       $AUGUSTUS_BIN_PATH,
                $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH
                );
            print LOG "Found script $string.\n" if ( $v > 3 );
            $perlCmdString = "perl $string $otherfilesDir/utr.gb.train "
                           . "$onlyTrainSize "
                           . "1> $otherfilesDir/randomSplit_utr2.log "
                           . "2> $errorfilesDir/randomSplit_utr2.err";
            print LOG "\n$perlCmdString\n" if ( $v > 3 );
            system("$perlCmdString") == 0 or die( "ERROR in file " . __FILE__
                . " at line " . __LINE__
                . "\nFailed to execute: $perlCmdString!\n" );
        }
        # changing UTR parameters in species config file to "on"
        print STDOUT "NEXT STEP: Setting value of \"UTR\" in "
            . "$AUGUSTUS_CONFIG_PATH/species/$species/$species\_parameters.cfg "
            . "to \"true\"\n";
        print LOG "\n\# "
            . (localtime)
            . ": Setting value of \"UTR\" in "
            . "$AUGUSTUS_CONFIG_PATH/species/$species/$species\_parameters.cfg "
            . "to \"true\"\n"
            if ( $v > 3 );
        setParInConfig(
            $AUGUSTUS_CONFIG_PATH
                . "/species/$species/$species\_parameters.cfg",
            "UTR", "on"
        );
        setParInConfig(
            $AUGUSTUS_CONFIG_PATH
                . "/species/$species/$species\_parameters.cfg",
            "print_utr", "on"
        );
    }


    if ( !uptodate( [ "utr.gb"], ["optimize.utr.out"] ) ) {

        # prepare metaparameter file
        my $metaUtrName = $species . "_metapars.utr.cfg";
        if (not( -e $AUGUSTUS_CONFIG_PATH . "/species/$species/$metaUtrName" )
            )
        {
            # copy from generic as template
            $cmdString = "cp $AUGUSTUS_CONFIG_PATH"
                       . "/species/generic/generic_metapars.utr.cfg $AUGUSTUS_CONFIG_PATH"
                       . "/species/$species/$metaUtrName";
            print LOG "Copying utr metaparameter template file:\n$cmdString\n"
                if ( $v > 3 );
            system("$cmdString") == 0 or die( "ERROR in file " . __FILE__
                . " at line " . __LINE__
                . "\nFailed to execute: $cmdString!\n" );
        }
        $string = find(
            "optimize_augustus.pl", $AUGUSTUS_BIN_PATH,
            $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH
        );
        print LOG "Found script $string.\n" if ( $v > 3 );
        if($onlyTrainSize == 0){
            $perlCmdString = "perl $string --rounds=$rounds --species=$species "
                           . "--trainOnlyUtr=1 "
                           . "--metapars=$AUGUSTUS_CONFIG_PATH"
                           . "/species/$species/$metaUtrName --cpus=$CPU "
                           . "$otherfilesDir/utr.gb.train "
                           . "--UTR=on > $otherfilesDir/optimize.utr.out";
        }else{
            $perlCmdString = "perl $string --rounds=$rounds --species=$species "
                           . "--trainOnlyUtr=1  "
                           . "--onlytrain=$otherfilesDir/utr.gb.train.train "
                           . "--metapars=$AUGUSTUS_CONFIG_PATH"
                           . "/species/$species/$metaUtrName --cpus=$CPU "
                           . "$otherfilesDir/utr.gb.train.test "
                           . "--UTR=on > $otherfilesDir/optimize.utr.out "
                           . "2> $errorfilesDir/optimize.utr.err"
        }
        print LOG "Now optimizing meta parameters of AUGUSTUS for the UTR "
            . "model:\n" if ( $v > 3 );
        print LOG "Running \"$perlCmdString\"..." if ( $v > 3 );
        system("$perlCmdString") == 0 or die( "ERROR in file " . __FILE__
            . " at line " . __LINE__
            . "\nFailed to execute: $perlCmdString!\n" );
    }
    else {
        print "Skipping UTR parameter optimization. Already up to date.\n";
    }
}

close(LOG);

