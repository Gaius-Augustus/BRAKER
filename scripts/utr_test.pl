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



train_utr()



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

  # store predictions without UTRs, revert later if UTR model does not improve
  # predictions
    print LOG "mv $otherfilesDir/augustus.gff "
        . "$otherfilesDir/augustus.noUtr.gff\n"
        if ( $v > 3 );
    move( "$otherfilesDir/augustus.gff",
        "$otherfilesDir/augustus.noUtr.gff" );
    print LOG
        "mv $otherfilesDir/augustus.gtf $otherfilesDir/augustus.noUtr.gtf\n"
        if ( $v > 3 );
    move( "$otherfilesDir/augustus.gtf",
        "$otherfilesDir/augustus.noUtr.gtf" );
    print LOG
        "mv $otherfilesDir/augustus.aa $otherfilesDir/augustus.noUtr.aa\n"
        if ( $v > 3 );
    move( "$otherfilesDir/augustus.aa", "$otherfilesDir/augustus.noUtr.aa" );

    # copy species parameter files, revert later if UTR model does not improve
    # predictions
    print LOG "\# "
        . (localtime)
        . ": Create backup of current species parameters:\n"
        if ( $v > 3 );
    for (
        (   "$species" . "_exon_probs.pbl",
            "$species" . "_igenic_probs.pbl",
            "$species" . "_intron_probs.pbl"
        )
        )
    {
        print LOG "cp $AUGUSTUS_CONFIG_PATH/species/$species/$_ "
            . "$AUGUSTUS_CONFIG_PATH/species/$species/$_.noUTR\n"
            if ( $v > 3 );
        copy(
            "$AUGUSTUS_CONFIG_PATH/species/$species/$_",
            "$AUGUSTUS_CONFIG_PATH/species/$species/$_.noUTR"
            )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCopy failed!\n" );
    }
    chdir($otherfilesDir)
        or die( "ERROR in file "
            . __FILE__
            . " at line "
            . __LINE__
            . "\nCould not change into directory $otherfilesDir!\n" );

    # search all start and stop codons from augustus.noUtr.gtf and write them
    # to the file stops.and.starts.gff
    if ( !uptodate( ["augustus.noUtr.gtf"], ["stops.and.starts.gff"] ) ) {
        print LOG "\# "
            . (localtime)
            . ": extracting all stop and start codons from augustus.noUtr.gtf "
            . "to stops.and.starts.gff\n"
            if ( $v > 3 );
        my %nonRedundantCodons;
        my @tmpGffLine;
        open( AUG, "<", "augustus.noUtr.gtf" )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not open file augustus.noUtr.gtf!\n" );
        while ( defined( my $i = <AUG> ) ) {

       # TODO: we are not dealing with redundancy, correctly. Discarding
       #       duplicates is not the optimal solution because later, we filter
       #       for genes where both codons have UTR models. However, at this
       #       point in time, we ignore this matter and hope that we are left
       #       with a sufficient number of training examples.
            if ( $i =~ /\t(start_codon|stop_codon)\t/ ) {
                @tmpGffLine = split( /\t/, $i );
                if (not(defined(
                            $nonRedundantCodons{
                                "$tmpGffLine[0]_$tmpGffLine[3]_$tmpGffLine[4]_$tmpGffLine[6]"
                            }
                        )
                    )
                    )
                {
                    $nonRedundantCodons{
                        "$tmpGffLine[0]_$tmpGffLine[3]_$tmpGffLine[4]_$tmpGffLine[6]"
                    } = $i;
                }
            }
        }
        close(AUG)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file augustus.noUtr.gtf!\n" );
        open( CODON, ">", "stops.and.starts.gff" )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not open file stops.and.starts.gff!\n" );
        foreach my $key ( keys %nonRedundantCodons ) {
            print CODON $nonRedundantCodons{$key};
        }
        close(CODON)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file stops.and.starts.gff!\n" );
    }

    if ( !uptodate( ["$hintsfile"], ["rnaseq.utr.hints"] ) ) {

      # TODO: currently, only using AT-AG, not AC-AG or any other splice site.
      #       Possibly extend to other splice patterns.
        print LOG "\# "
            . (localtime)
            . ": filtering RNA-Seq hints for valid splice site AT-AG, storing "
            . "in rnsaeq.utr.hints\n"
            if ( $v > 3 );
        my %genome_hash;
        my $hash_key;
        open( FASTA, "<", $genome )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not open file $genome!\n" );
    LINE: while ( my $line = <FASTA> ) {
            next LINE if $line =~ m/^#/;    #discard comments
            if ( $line =~ /^>/ ) {
                chomp($line);
                $hash_key = substr( $line, 1, length($line) - 1 );
            }
            else {
                $line =~ s/[\x0A\x0D]+//g;
                $line =~ s/(\s+)(\n)(\r)//g;
                $line = uc($line);
                $genome_hash{$hash_key} .= $line;
            }
        }
        close(FASTA)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file $genome!\n" );
        open( HINTS, "<", $hintsfile )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not open file $hintsfile!\n" );
        my @gff;
        my ( $siteA, $siteB, $given, $lastCol );
        my $splice = "ATAG";
        open( UTRHINTS, ">", "rnaseq.utr.hints" )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not open file rnaseq.utr.hints!\n" );
    LINE: while ( my $line = <HINTS> ) {
            @gff = split( /\t/, $line );
            if ( ( $gff[1] eq "b2h" ) && ( $gff[2] eq "intron" ) )
            {    # make sure to use only intron hints from RNA-Seq data
                $siteA
                    = substr( $genome_hash{ $gff[0] }, ( $gff[3] - 1 ), 2 );
                $siteB
                    = substr( $genome_hash{ $gff[0] }, ( $gff[4] - 2 ), 2 );
                $given = $siteA . $siteB;
                if ( $gff[8] =~ m/mult=(\d+)/ ) {
                    $lastCol = "mult=$1_$splice\n";
                }
                else {
                    $lastCol = "mult=1_$splice\n";
                }
                if ( uc($given) =~ m/$splice/ ) {
                    print $gff[0] . "\t"
                        . $gff[1] . "\t"
                        . $gff[2] . "\t"
                        . $gff[3] . "\t"
                        . $gff[4] . "\t"
                        . $gff[5] . "\t+\t"
                        . $gff[7] . "\t"
                        . $lastCol;
                }
                else {
                    $given = reverse $given;
                    $given =~ tr/ACGTacgt/TGCAtgca/;
                    if ( uc($given) =~ m/$splice/ ) {
                        print $gff[0] . "\t"
                            . $gff[1] . "\t"
                            . $gff[2] . "\t"
                            . $gff[3] . "\t"
                            . $gff[4] . "\t"
                            . $gff[5] . "\t-\t"
                            . $gff[7] . "\t"
                            . $lastCol;
                    }
                }
            }
        }
        close(UTRHINTS)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file rnaseq.utr.hints!\n" );
        close(HINTS)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file $hintsfile!\n" );
    }

    # create wiggle file from bam files
    if ( !uptodate( ["$hintsfile"], ["rnaseq.wig"] ) ) {
        if ( scalar(@bam) > 1 ) {
            print LOG "\# "
                . (localtime)
                . ": converting bam files to wiggle file rnaseq.wig\n"
                if ( $v > 3 );
            $cmdString = "";
            if ($nice) {
                $cmdString .= "nice ";
            }
            $cmdString .= "$BAMTOOLS_BIN_PATH/bamtools merge ";
            foreach (@bam) {
                chomp;
                $cmdString .= "-in $_ ";
            }
            $cmdString .= "-out merged.bam";
            print LOG "\n$cmdString\n\n" if ( $v > 3 );
            system("$cmdString")
                or die( "ERROR in file "
                    . __FILE__
                    . " at line "
                    . __LINE__
                    . "\nFailed to execute: $cmdString!\n" );
        }
        else {
            print LOG "\# "
                . (localtime)
                . ":  Creating softlink to bam file $bam[0]...\n"
                if ( $v > 3 );
            $cmdString = "ln -s $bam[0] merged.bam";
            print LOG "$cmdString\n" if ( $v > 3 );
            system($cmdString) == 0
                or die( "ERROR in file "
                    . __FILE__
                    . " at line "
                    . __LINE__
                    . "\nFailed to execute: $cmdString!\n" );
        }
        print LOG "\# " . (localtime) . ": Creating wiggle file...\n"
            if ( $v > 3 );
        $cmdString = "";
        if ($nice) {
            $cmdString .= "nice ";
        }
        $cmdString
            .= "$bam2wigPath merged.bam >merged.wig 2> $otherfilesDir/bam2wig.err";
        print LOG "\n$cmdString\n" if ( $v > 3 );
        system("$cmdString") == 0
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nFailed to execute: $cmdString!\n" );
    }

    # call utrrnaseq
    if ( !uptodate( ["$hintsfile"], ["utrs.gff"] ) ) {
        print LOG "\# " . (localtime) . ": Creating utrs.gff\n" if ( $v > 3 );
        $cmdString = "";
        if ($nice) {
            $cmdString .= "nice ";
        }
        $cmdString
            .= "$rnaseq2utrPath -G $genome -O stops.and.starts.gff -I rnaseq.utr.hints -W rnaseq.wig --outFileName=utrs.gff $rnaseq2utr_args 2> $otherfilesDir/rnaseq2utr.err";
        print LOG "\n$cmdString\n" if ( $v > 3 );
        system("$cmdString") == 0
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nFailed to execute: $cmdString!\n" );
    }

    # create genbank file with genes that have to utrs
    if (!uptodate(
            [ "utrs.gff",    "augustus.noUtr.gtf" ],
            [ "bothutr.lst", "bothutr.test.gb" ]
        )
        )
    {
        print LOG "\# "
            . (localtime)
            . ": Creating gb file for UTR training\n"
            if ( $v > 3 );

        # extract subset of genes, where we have both UTRs
        open( UTR, "<", "utrs.gff" )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCan not open file utrs.gff!\n" );
        open( TRLST, ">", "tr.lst" );
        while (<UTR>) {
            s/.*\t(\S+UTR)\t.*transcript_id \"(\S+)\".*/$2\t$1/;
            print TRLST;
        }
        close(UTR)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file utrs.gff!\n" );
        close(TRLST)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file tr.lst!\n" );
        $cmdString = "";
        if ($nice) {
            $cmdString .= "nice ";
        }
        $cmdString .= "cat tr.lst | sort -u > tr_temp.lst";
        print LOG "\n$cmdString\n" if ( $v > 3 );
        system($cmdString) == 0
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nFailed not execute $cmdString!\n" );
        print LOG "\nrm tr.lst\n" if ( $v > 3 );
        unlink("tr.lst");
        $cmdString = "mv tr_temp.lst tr.lst";
        print LOG "\n$cmdString\n" if ( $v > 3 );
        system($cmdString) == 0
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nFailed not execute $cmdString!\n" );
        open( TR, "tr.lst" )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCan not open file tr.lst!\n" );
        open( BOTH, ">", "bothutr.lst" );
        my $Fld1;
        my $prev;

        while (<TR>) {
            ($Fld1) = split( '\t', $_, -1 );
            if ( $Fld1 eq $prev ) {
                print BOTH "$prev\n";
            }
        }
        close(TR)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file tr.lst!\n" );
        close(BOTH)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file bothutr.lst!\n" );
        $cmdString = "";
        if ($nice) {
            $cmdString .= "nice ";
        }
        $cmdString .= "cat utrs.gff augustus.noUtr.gtf > genes.gtf_temp";
        print LOG "\n$cmdString\n" if ( $v > 3 );
        system($cmdString) == 0
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nFailed not execute $cmdString!\n" );
        open( GENES, "<", "genes.gtf_temp" )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCan not open the file genes.gtf_temp!\n" );
        open( WRITEGENES, ">", "genes.gtf_unsort" );
        while (<GENES>) {
            if (/(CDS|UTR)\t/) {
                print WRITEGENES "Not sure what\n";
            }
        }
        close(GENES)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file genes.gtf_temp!\n" );
        close(WRITEGENES)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file genes.gtf_unsort!\n" );
        $cmdString = "";
        if ($nice) {
            $cmdString .= "nice ";
        }
        $cmdString
            .= "cat genes.gtf_unsort | sort -n -k 4,4 | sort -s -k 10,10 | sort -s -k 1,1 > genes.gtf";
        print LOG "\n$cmdString\n" if ( $v > 3 );
        system($cmdString) == 0
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nFailed not execute $cmdString!\n" );
        $string = find(
            "gff2gbSmallDNA.pl",    $AUGUSTUS_BIN_PATH,
            $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH
        );
        $perlCmdString = "";
        if ($nice) {
            $perlCmdString .= "nice ";
        }
        $perlCmdString
            .= "perl $string genes.gtf $genome $flanking_DNA bothutr.test.gb --good=bothutr.lst 1> $otherfilesDir/gff2gbSmallDNA.utr.stdout "
            . "2> $otherfilesDir/gff2gbSmallDNA.utr.stderr";
        print LOG "\n$perlCmdString\n" if ( $v > 3 );
        system("$perlCmdString") == 0
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nFailed to execute: $perlCmdString!\n" );
    }

# create train.gb and onlytrain.gb
# train.gb contains a small proportion of genes with utr and up to 150 genes with UTR
# onlytrain.gb contains all other genes (with and without utr)
    if (!uptodate(
            [ "bothutr.test.gb", "../training.gb.train.test" ],
            [ "train.gb",        "onlytrain.gb" ]
        )
        )
    {
        # evaluate m: size of smaller set of utr examples
        my $m;

        # count the block number in bothutr.test.gb
        my $count = 0;
        open( TEMP1, "<", "bothutr.test.gb" )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCan not open the file bothutr.test.gb! \n" );
        while (<TEMP1>) {
            $count++ if ( $_ =~ m/LOCUS/ );
        }
        close(TEMP1)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file bothutr.test.gb!\n" );
        if   ( $count >= 150 ) { $m = 150 }
        else                   { $m = $count }
        if ( $count < 50 ) {
            die(      "ERROR in file "
                    . __FILE__
                    . " at line "
                    . __LINE__
                    . "\n Number of UTR training examples is smaller than 50. Abort "
                    . "UTR training. If this is the only error message, the "
                    . "AUGUSTUS parameters for your species were optimized ok, "
                    . "but you are lacking UTR parameters. Do not attempt to "
                    . "predict genes with UTRs for this species using the current "
                    . "parameter set!\n" );
            exit(1);
        }

        # evaluate n: size of smaller set of no utr examples
        my $n;

        # count the block number in training.gb.train.test
        $count = 0;
        open( TEMP2, "train.gb.test" )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCan not open the file train.gb.test!\n" );
        while (<TEMP2>) {
            $count++ if ( $_ =~ m/LOCUS/ );
        }
        close(TEMP2)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file train.gb.test!\n" );
        if ( $count >= 50 ) {
            $n = 50;
        }
        else {
            $n = $count;
        }

        # extract traininging set for UTR model
        $string = find(
            "randomSplit.pl",       $AUGUSTUS_BIN_PATH,
            $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH
        );
        print LOG "Found script $string.\n" if ( $v > 3 );
        $perlCmdString = "perl $string bothutr.test.gb $m";
        print LOG "\nperlCmdString\n" if ( $v > 3 );
        system("$perlCmdString") == 0
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nFailed to execute: $perlCmdString!\n" );
        $perlCmdString = "perl $string train.gb.test $n";
        print LOG "\n$perlCmdString\n" if ( $v > 3 );
        system("$perlCmdString") == 0
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nFailed to execute: $perlCmdString!\n" );
        my $delete;
        open( GB, "<", "train.gb.test.test" )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCan not open file train.gb.test.test!\n" );
        open( NOMRNA, ">", "nomrna.test.gb" )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not open file nomrna.test.gb!\n" );

        while (<GB>) {
            $delete = 1 if /mRNA/;
            $delete = 0 if /CDS/;
            print NOMRNA if ( !$delete );
        }
        close(GB)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file train.gb.test.test!\n" );
        close(NOMRNA)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file nomrna.test.gb!\n" );
        $cmdString = "cat nomrna.test.gb bothutr.test.gb.test > train.gb";
        print LOG "\n$cmdString\n" if ( $v > 3 );
        system("$cmdString") == 0
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nFailed to execute: $cmdString\n" );

        # count how many genes are contained in train.gb
        my $counter_gen = 0;
        open( TS, "train.gb" );
        while (<TS>) {
            $counter_gen++ if (/^     CDS             /);
        }
        close(TS)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file train.gb!\n" );
        print LOG "Have constructed a training set train.gb for UTRs with "
            . "$counter_gen genes\n"
            if ( $v > 3 );
        print LOG "Deleting nomrna.test.gb, train.gb.test.test, "
            . "train.gb.test.train\n"
            if ( $v > 3 );
        unlink("nomrna.test.gb");
        unlink("train.gb.test.test");
        unlink("train.gb.test.train");

        # create onlytrain training set only used for training #
        open( ONLYTRAIN, "<", "train.gb.train" )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not open file train.gb.train!\n" );
        open( CDSONLY, ">", "cdsonly.gb" )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not open file cdsonly.gb!\n" );

        # delete the mRNA part up to the next CDS tag
        $delete = 0;
        while (<ONLYTRAIN>) {
            $delete = 1 if /mRNA/;
            $delete = 0 if /CDS/;
            print CDSONLY if ( !$delete );
        }
        close(ONLYTRAIN)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file train.gb.train!\n" );
        close(CDSONLY)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file cdsonlyl.gb!\n" );

       # construct the disjoint sets: remove training UTR genes from onlytrain
       # UTR gene set (train.utronly.gb)
        open( TRAIN, "<", "train.gb" )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not open the file train.gb!\n" );
        open( REMOVE, ">", "remove.lst" )
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not open file remove.lst!\n" );
        my $locustag = 0;
        while (<TRAIN>) {
            if (m/LOCUS\s+(\S+)_\d+-\d+/) {
                $locustag = 0;
                print REMOVE "$1_";
            }
            elsif (m/gene="(\S+)\.t\d+/) {
                if ( $locustag == 0 ) {
                    print REMOVE $1 . "\n";
                }
                $locustag = 1;
            }
        }
        close(TRAIN)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file train.gb!\n" );
        close(REMOVE)
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nCould not close file remove.lst!\n" );
        $string = find(
            "filterGenes.pl",       $AUGUSTUS_BIN_PATH,
            $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH
        );
        print LOG "Found script $string.\n" if ( $v > 3 );
        $perlCmdString
            = "perl $string remove.lst bothutr.test.gb > train.utronly.gb";
        print LOG "\n$perlCmdString\n" if ( $v > 3 );
        system("$perlCmdString") == 0
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nFailed to execute: $perlCmdString!\n" );
        $cmdString = "cat cdsonly.gb train.utronly.gb > onlytrain.gb";
        print LOG "\n$cmdString\n" if ( $v > 3 );
        system("$cmdString") == 0
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nfailed to execute: $cmdString!\n" );

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
    if ( !uptodate( [ "train.gb", "onlytrain.gb" ], ["optimize.utr.out"] ) ) {

        # prepare metaparameter file
        my $metaUtrName = $species . "_metapars.utr.cfg";
        if (not( -e $AUGUSTUS_CONFIG_PATH . "/species/$species/$metaUtrName" )
            )
        {
            # copy from generic as template
            $cmdString
                = "cp $AUGUSTUS_CONFIG_PATH"
                . "/species/generic/generic_metapars.utr.cfg $AUGUSTUS_CONFIG_PATH"
                . "/species/$species/$metaUtrName";
            print LOG
                "Copying utr metaparameter template file:\n$cmdString\n"
                if ( $v > 3 );
            system("$cmdString") == 0
                or die( "ERROR in file "
                    . __FILE__
                    . " at line "
                    . __LINE__
                    . "\nFailed to execute: $cmdString!\n" );
        }
        $string = find(
            "optimize_augustus.pl", $AUGUSTUS_BIN_PATH,
            $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH
        );
        print LOG "Found script $string.\n" if ( $v > 3 );
        $perlCmdString
            = "perl $string --rounds=$rounds --species=$species --trainOnlyUtr=1 "
            . "--onlytrain=onlytrain.gb  --metapars=$AUGUSTUS_CONFIG_PATH"
            . "/species/$species/$metaUtrName train.gb --UTR=on > optimize.utr.out";
        print LOG
            "Now optimizing meta parameters of AUGUSTUS for the UTR model:\n"
            if ( $v > 3 );
        print LOG "Running \"$perlCmdString\"..." if ( $v > 3 );
        system("$perlCmdString") == 0
            or die( "ERROR in file "
                . __FILE__
                . " at line "
                . __LINE__
                . "\nFailed to execute: $perlCmdString!\n" );
    }
    else {
        print "Skipping UTR parameter optimization. Already up to date.\n";
    }
}
