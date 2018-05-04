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
my $flanking_DNA;
my $rounds = 3;
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
                       .  "&> $errorfilesDir/bam.merge.err";
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
                   .  "$rnaseq2utr_args 2> $errorfilesDir/rnaseq2utr.err";
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
        open( TRLST, ">", "$otherfilesDir/tr.lst" ) or die ("ERROR in file "
            . __FILE__ . " at line ". __LINE__
            . "\nCan not open file $otherfilesDir/tr.lst!\n");
        while (<UTR>) {
            $_ = ~ s/.*\t(\S+UTR)\t.*transcript_id \"(\S+)\".*/$2\t$1/;
            print TRLST $_;
        }
        close(UTR) or die( "ERROR in file " . __FILE__ . " at line " . __LINE__
                . "\nCould not close file $otherfilesDir/utrs.gff!\n" );
        close(TRLST) or die( "ERROR in file " . __FILE__ . " at line "
            . __LINE__ . "\nCould not close file $otherfilesDir/tr.lst!\n" );
        $cmdString = "";
        if ($nice) {
            $cmdString .= "nice ";
        }
        $cmdString .= "cat $otherfilesDir/tr.lst | sort -u > "
                   .  "$otherfilesDir/tr_temp.lst "
                   .  "2> $errorfilesDir/tr_tmp_sort.err";
        print LOG "\n$cmdString\n" if ( $v > 3 );
        system($cmdString) == 0 or die( "ERROR in file " . __FILE__
            . " at line " . __LINE__ . "\nFailed not execute $cmdString!\n" );
        print LOG "\nrm $otherfilesDir/tr.lst\n" if ( $v > 3 );
        unlink("$otherfilesDir/tr.lst");
        $cmdString = "mv $otherfilesDir/tr_temp.lst $otherfilesDir/tr.lst";
        print LOG "\n$cmdString\n" if ( $v > 3 );
        system($cmdString) == 0 or die( "ERROR in file " . __FILE__
            . " at line " . __LINE__ . "\nFailed not execute $cmdString!\n" );

        open( TR, "<", "$otherfilesDir/tr.lst" ) or die( "ERROR in file "
            . __FILE__ . " at line " . __LINE__
            . "\nCan not open file $otherfilesDir/tr.lst!\n" );
        my %both;
        while (<TR>) {
            my @t = split(/\t/);
            $t[8]=~m/transcript_id \"\S+\"";
            my $txid = $1;
            if($t[2] =~ m/UTR/){
                $both{$txid}{$t[2]} = $_;
            }
        }
        close(TR) or die( "ERROR in file " . __FILE__ . " at line " . __LINE__
            . "\nCould not close file $otherfilesDir/tr.lst!\n" );

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
                   .  "$otherfilesDir/augustus.hints.noUtr.gtf "
                   .  "> $otherfilesDir/genes.gtf_temp "
                   .  "2> $errorfilesDir/cat_utrs_augustus_noUtrs.err";
        print LOG "\n$cmdString\n" if ( $v > 3 );
        system($cmdString) == 0 or die( "ERROR in file " . __FILE__
            . " at line " . __LINE__
            . "\nFailed not execute $cmdString!\n" );

        open( GENES, "<", "$otherfilesDir/genes.gtf_temp" )
            or die( "ERROR in file " . __FILE__ . " at line " . __LINE__
                . "\nCan not open the file $otherfilesDir/genes.gtf_temp!\n" );
        open( WRITEGENES, ">", "$otherfilesDir/genes.gtf_unsort" );
        while (<GENES>) {
            if (/(CDS|UTR)\t/) {
                print WRITEGENES "Not sure what\n"; # TODO: FIX!
            }
        }
        close(GENES) or die( "ERROR in file " . __FILE__ . " at line "
            . __LINE__ . "\nCould not close file "
            . "$otherfilesDir/genes.gtf_temp!\n" );
        close(WRITEGENES) or die( "ERROR in file " . __FILE__ . " at line "
            . __LINE__
            . "\nCould not close file $otherfilesDir/genes.gtf_unsort!\n" );
        $cmdString = "";
        if ($nice) {
            $cmdString .= "nice ";
        }
        $cmdString .= "cat genes.gtf_unsort | sort -n -k 4,4 | "
                   .  "sort -s -k 10,10 | sort -s -k 1,1 > genes.gtf "
                   .  "2> $errorfilesDir/genes.gtf_unsort.err";
        print LOG "\n$cmdString\n" if ( $v > 3 );
        system($cmdString) == 0 or die( "ERROR in file " . __FILE__
            . " at line " . __LINE__ . "\nFailed not execute $cmdString!\n" );

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
                       . "&> $errorfilesDir/randomSplit_utr1.err";
        print LOG "\nperlCmdString\n" if ( $v > 3 );
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
                           . "&> $errorfilesDir/randomSplit_utr2.err";
            print LOG "\nperlCmdString\n" if ( $v > 3 );
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
                           . "--trainOnlyUtr=1  "
                           . "--metapars=$AUGUSTUS_CONFIG_PATH"
                           . "/species/$species/$metaUtrName "
                           . "$otherfilesDir/utr.gb.train "
                           . "--UTR=on > $otherfilesDir/optimize.utr.out";
        }else{
            $perlCmdString = "perl $string --rounds=$rounds --species=$species "
                           . "--trainOnlyUtr=1  "
                           . "--onlytrain=$otherfilesDir/utr.gb.train.train"
                           . "--metapars=$AUGUSTUS_CONFIG_PATH"
                           . "/species/$species/$metaUtrName "
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

