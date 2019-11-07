#!/usr/bin/env perl

####################################################################################################
#                                                                                                  #
# startAlign.pl - run alignment program Spaln2, Exonerate or GenomeThreader to align protein       #
#                  sequences to genome                                                             #
#                                                                                                  #
# Author: Katharina Hoff & Simone Lange                                                            #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# Last changes: Feb 15th 2018                                                                      #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
#                                                                                                  #
####################################################################################################

# ----------------------------------------------------------------
# | clean up after merging          | Katharina Hoff |15.02.2018 |
# | added path to alignment tools   | Katharina Hoff |12.07.2017 |
# | changed script description      | Katharina Hoff |11.07.2017 |
# | first outline                   | Simone Lange   |28.10.2015 |
# | changed default from region to  |                |01.12.2015 |
# | whole sequence (noreg to reg)   |                |           |
# | added "wholeProt option" (just  |                |           |
# | leave out list and pos files)   |                |           |
# ----------------------------------------------------------------

use Getopt::Long;
use Cwd;
use File::Path qw(make_path rmtree);
use File::Spec::Functions qw(rel2abs);
use Parallel::ForkManager;

use strict;
use warnings;

my $usage = <<'ENDUSAGE';

startAlign.pl split genome file in single sequences or sequence parts and protein file according
to contigIDs then run alignment program exonerate, spaln or genomeThreader (gth)
for each contigID or each sequence. When no list and/or pos file(s) are/is assigned
the program will use the whole protein file.

SYNOPSIS

startAlign.pl [OPTIONS] --genome=genome.fa --prot=db.fa  --prg=gth|exonerate|spaln

--genome=genome.fa          fasta file with DNA sequences
--prot=db.fa                fasta file with protein sequences


OPTIONS

--help                       Print this help message
--CPU=n                      Specifies the maximum number of CPUs that can be used during
WARNING: using more than one CPU will produce harddrive read/write
processes that may be speed limiting in case of GenomeThreader
--dir=path/to/dir            Set path to working directory. In the working directory results
and temporary files are stored.
--list=BLAST.hit.list        Contains contig and protein ID. Format: contigID proteinID
--log=startAlign.log         Log file
--maxintronlen=n             Exonerate option: Alignments with longer gaps are discarded (default 30000).
--reg                        Use region parts and not whole sequences.
--offset=n                   This many bp are added before and after cutout coordinates.
--prg=s                      Alignment program to call. Valid options are 'gth', 'spaln' or 'exonerate'.
--pos=dna.pos                Contains information on contigs and genome sequence. Format
contigID nr_of_prots_mapped start end strand chrID
--ALIGNMENT_TOOL_PATH=/path/to/binary
Path to alignment tool binary, either exonerate or Splan or Genome Threader.
By default, if no path is given, script assumes they are in the current
$PATH bash variable.
--args=s                     additional command line parameters for alignment tool to be executed,
                             example: --args="-prinmatchlen 24 -prseedlength 10 -prhdist 4"
--nice                       Execute all system calls within braker.pl and its submodules with bash "nice"
(default nice value)



DESCRIPTION

Example:

startAlign.pl [OPTIONS] --genome=genome.fa --prot=db.fa --list=BLAST.hit.list --pos=dna.pos --prg=gth

ENDUSAGE

my $alnArgs;
my $alignDir;    # directory for alignment output files
my $errorfile;   # error file name
my $cmdString;   # to store shell commands
my %contigIDs;   # hash for contig IDs
my $counterW
    = 0;         # print message 'add amino ...' for gth only once per protein
my $CPU = 1;     # number of CPUs that can be used
my $dir;         # working superdirectory where program is called from
my $genome_file; # genome file
my $help;        # print usage
my $list_file;   # blast hit list file
my $log;         # log file name
my $maxintronlen = 30000;   # maximal intron length (for usage with exonerate)
my $offset       = 10000;   # offset for cutout [start-offset, end+offset]
my $pos_file;          # position list file
my $prgsrc;            # source program (exonerate, spaln or gth)
my $prot_file;         # protein database file
my $prot_file_base;    # protein database file name base for $protwhole option
my $prot_addstop_file
    ;    # file name that will be automatically derived from $prot_file_base
my %protIDs;    # hash for protein IDs
my $protWhole = 0;    # use whole protein database file
my $reg       = 0;    # use regions
my %seq;              # hash for genome sequences
my $stdoutfile;       # standard output file name
my $tmpDir;    # temporary directory for storing protein and genome part files
my $nice; # flag that determines whether system calls should be executed with bash nice
          #(default nice value)

my $spalnErrAdj = 80; # version spaln2.2.0 misses a line while "counting coordinates" (+80 because of how the fasta files are printed here, see line 529-534); error still persists with version 2.3.1

# gth options           # and other values that have proven to work well for Drosophila on chromosome 2L
my $gcmincoverage = 80;    # 70 - 95
my $prhdist       = 2;     # 2 - 6
my $prseedlength = 20;   # 19, 20, 21 (gth: error: -prminmatchlen must be >= -prseedlength)
my $prminmatchlen = 20;
my $ALIGNMENT_TOOL_PATH;

if ( @ARGV == 0 ) {
    print "$usage\n";
    exit(0);
}

GetOptions(
    'CPU=i'                 => \$CPU,
    'dir=s'                 => \$dir,
    'genome=s'              => \$genome_file,
    'list=s'                => \$list_file,
    'log=s'                 => \$log,
    'maxintron=i'           => \$maxintronlen,
    'reg!'                  => \$reg,
    'offset=i'              => \$offset,
    'pos=s'                 => \$pos_file,
    'prg=s'                 => \$prgsrc,
    'prot=s'                => \$prot_file,
    'help!'                 => \$help,
    'nice!'                 => \$nice,
    'ALIGNMENT_TOOL_PATH=s' => \$ALIGNMENT_TOOL_PATH,
    'args=s'                => \$alnArgs
    );

if ($help) {
    print $usage;
    exit(0);
}

# if no working directory is set, use current directory
if ( !defined($dir) ) {
    $dir = cwd();
}

# remove tailing slash from ALIGNMENT_TOOL_PATH
if ( defined($ALIGNMENT_TOOL_PATH) ) {
    my $last_char = substr( $ALIGNMENT_TOOL_PATH, -1 );
    if ( $last_char eq "\/" ) {
        chop($ALIGNMENT_TOOL_PATH);
    }
}

# check whether position and list files are specified
if ( defined($pos_file) && defined($list_file) ) {
    # check whether position file exists
    if ( !-e $pos_file ) {
        print STDERR
        "ERROR: in file " . __FILE__ ." at line ". __LINE__ ."\nPosition file $pos_file does not exist. Please check.\n";
        exit(1);
    } else {
        $pos_file = rel2abs($pos_file);
    }

    # check whether list file exists
    if ( !-e $list_file ) {
        print STDERR
        "ERROR: in file " . __FILE__ ." at line ". __LINE__ ."\nList file $list_file does not exist. Please check.\n";
        exit(1);
    } else {
        $list_file = rel2abs($list_file);
    }
} else {
    $protWhole = 1;
}

if ( !defined($prgsrc) ) {
    print STDERR
    "ERROR: in file " . __FILE__ ." at line ". __LINE__ ."\nPlease assign the source program with --prg. Possible Options are 'exonerate', 'spaln' or 'gth'.\n";
    exit(1);
}

# check program source option
if ( $prgsrc ne "exonerate" && $prgsrc ne "spaln" && $prgsrc ne "gth" ) {
    print STDERR
    "ERROR: in file " . __FILE__ ." at line ". __LINE__ ."\nInvalid value '$prgsrc' for option --prg. Possible Options are 'exonerate', 'spaln' or 'gth'. Please check.\n";
    exit(1);
}

if ( !defined($log) ) {
    $log = "$dir/startAlign_$prgsrc.log";
}
open( LOG, ">" . $log ) or die "ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCannot open file: $log\n";

$tmpDir   = "$dir/tmp_$prgsrc";
$alignDir = "$dir/align_$prgsrc";
if ( ! (-d $tmpDir) ) {
    make_path($tmpDir) or die ("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCannot create directory $tmpDir!\n");
    print LOG "\# " . (localtime) . ": create working directory $tmpDir\n";
    print LOG "mkdir $tmpDir\n\n";
} else {
    $cmdString = "rm -r $tmpDir";
    print LOG "\# " . (localtime) . ": delete existing files from $tmpDir\n";
    print LOG "$cmdString\n\n";
    system("$cmdString") == 0 or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nfailed to execute: $cmdString!\n");
    make_path($tmpDir);
}

# check whether genome file is specified
if ( defined($genome_file) ) {

    # check whether genome file exists
    if ( !-e $genome_file ) {
        print STDERR "ERROR: in file " . __FILE__ ." at line ". __LINE__ ."\nGenome file $genome_file does not exist!\n";
        exit(1);
    } else {
        $genome_file = rel2abs($genome_file);
        my $nDots = ( $genome_file =~ tr/\.// );
        my $linkName = $genome_file;
        if ( $nDots > 1 ) {    # more than one dot in file name, need to avoid that for spaln
            for ( my $i = 0; $i < ( $nDots - 1 ); $i++ ) {
                $linkName =~ s/\./_/;
            }
        }

        # create link to genome file in tmpDir (Spaln is looking for related database files, there)
        my @pathParts = split( /\//, $linkName );
        $linkName = $tmpDir . "/" . $pathParts[ scalar(@pathParts) - 1 ];
        if ( !-e $linkName ) {
            $cmdString = "ln -s $genome_file $linkName";
            system($cmdString) == 0 or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nFailed to execute $cmdString!\n");
            print LOG "\# " . (localtime) . ": $cmdString\n";
        } else {
            print LOG "\# "
                . (localtime)
                . " WARNING: $linkName does already exist. Will assume that $linkName is a link to $genome_file!\n";
        }
        $genome_file = $linkName;
    }
} else {
    print STDERR "ERROR: in file " . __FILE__ ." at line ". __LINE__ ."\nNo genome file specified!\n";
    exit(1);
}

# check whether protein file is specified
if ( defined($prot_file) ) {

    # check whether protein file exists
    if ( !-e $prot_file ) {
        print STDERR "ERROR: in file " . __FILE__ ." at line ". __LINE__ ."\nProtein file $prot_file does not exist!\n";
        exit(1);
    } else {
        $prot_file = rel2abs($prot_file);
        my $linkName = $prot_file;
        my $nDots = ( $prot_file =~ tr/\.// );
        if ( $nDots > 1 ) {    # more than one dot in file name, need to avoid that for spaln
            for ( my $i = 0; $i < ( $nDots - 1 ); $i++ ) {
                $linkName =~ s/\./_/;
            }
        }
        my @pathParts = split( /\//, $linkName );
        $linkName = $tmpDir . "/" . $pathParts[ scalar(@pathParts) - 1 ];
        if ( !-e $linkName ) {
            $cmdString = "ln -s $prot_file $linkName";
            system($cmdString) == 0 or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nFailed to execute $cmdString!\n");
            print LOG "\# " . (localtime) . ": $cmdString\n";
        } else {
            print LOG "\# "
                . (localtime)
                . " WARNING: $linkName does already exist. Will assume that $linkName is a link to $prot_file!\n";
        }
        $prot_file = $linkName;
        my @a = split( /\//, $prot_file );
        $prot_file_base = $a[-1];

        # define prot_file.addstop file name
        $prot_addstop_file = $prot_file_base;
        @pathParts = split( /\./, $prot_file_base );
        if ( scalar(@pathParts) == 2 ) {
            $prot_addstop_file = $pathParts[0] . "_addstop." . $pathParts[1];
        } else {
            die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nUnexpected error in prot_addstop_file name definition\n");
        }
    }
} else {
    print STDERR "ERROR: in file " . __FILE__ ." at line ". __LINE__ ."\nNo protein file specified!\n";
    exit(1);
}

my $last_char = substr( $dir, -1 );
if ( $last_char eq "\/" ) {
    chop($dir);
}

if ( $prgsrc eq "spaln" ) {
    # check for spaln2 environment variables
    if ( !$ENV{'ALN_DBS'} ) {
        print STDERR
        "ERROR: in file " . __FILE__ ." at line ". __LINE__ ."\nThe environment variable ALN_DBS for spaln2 is not defined. Please export an environment variable with:' export ALN_DBS=/path/to/spaln2/seqdb'\n";
        exit(1);
    }
    if ( !$ENV{'ALN_TAB'} ) {
        print STDERR
        "ERROR: in file " . __FILE__ ." at line ". __LINE__ ."\nThe environment variable ALN_TAB for spaln2 is not defined. Please export an environment variable with:' export ALN_TAB=/path/to/spaln2/table'\n";
        exit(1);
    }
}

# add 80 bp to coordinates, if --prg is spaln, 0 otherwise
if ( $prgsrc ne "spaln" ) {
    $spalnErrAdj = 0;
}

# read in list and
if ( !$protWhole ) {
    read_files();
}

# get protein sequences from files and store in hash
prots();

if ( $CPU > 1 || ( $CPU == 1 && !$protWhole ) ) {
    get_seqs();
}

start_align();
clean_up();

close(LOG) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close file $log!\n");

############### sub functions ##############

# get protein sequences from files
sub read_files {
    open( POS, $pos_file ) or die "ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCannot open file: $pos_file\n";
    print LOG "\# " . (localtime) . ": read in positions from $pos_file\n";
    while (<POS>) {
        chomp;
        my @line = split( /[\s, ]/, $_ );
        if ( scalar(@line) == 6 ) {
            $contigIDs{ $line[0] }{"start"} = $line[2] - 1;
            $contigIDs{ $line[0] }{"end"}   = $line[3] - 1;
            my @seqname = split( /[\s, ]/, $line[5] );    # seqname only up to first white space
            $contigIDs{ $line[0] }{"seq"} = $seqname[0];
        }
    }
    close(POS) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close file $pos_file!\n");

    open( LIST, $list_file ) or die "ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCannot open file: $list_file\n";
    print LOG "\# " . (localtime) . ": read in positions from $list_file\n";
    while (<LIST>) {
        chomp;
        my @line = split( /[\s+, ]/, $_ );
        if ( scalar(@line) == 2 ) {
            if ($reg == 1) {
                push( @{ $protIDs{ $line[1] } }, $line[0] );
            } else {
                if ( defined( $contigIDs{ $line[0] }{"seq"} ) ) {
                    push(@{ $protIDs{ $line[1] } }, $contigIDs{ $line[0] }{"seq"}) unless grep { $_ eq $contigIDs{ $line[0] }{"seq"} } @{ $protIDs{ $line[1] } };
                }
            }
        }
    }
    close(LIST) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close file $list_file!\n");
}

sub get_seqs {
    open( FASTA, $genome_file ) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCannot open file: $genome_file\n");
    print LOG "\# "
        . (localtime)
        . ": read in DNA sequence from $genome_file\n";
    $/ = ">";
    while (<FASTA>) {
        s/>$//;    # see getAnnoFasta.pl
        next unless m/\S+/;    # see getAnnoFasta.pl
        /(.*)\n/;              # see getAnnoFasta.pl
        my $seqname      = $1; # see getAnnoFasta.pl
        my $sequencepart = $'; # see getAnnoFasta.pl
        $seqname =~ s/\s.*//;    # seqname only up to first white space (see getAnnoFasta.pl)
        $sequencepart =~ s/\s//g;    # see getAnnoFasta.pl
        $seq{$seqname} = $sequencepart;
    }
    $/ = "\n";
    close(FASTA) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close file $genome_file!\n");
}

sub prots {
    open( PROT, $prot_file ) or die "ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCannot open file: $prot_file\n";
    print LOG "\# "
        . (localtime)
        . ": read in fasta sequences from $prot_file\n";
    my $seqname;
    my $sequencepart = "";
    while (<PROT>) {
        if ( $_ =~ m/^>/ ) {
            chomp;
            if ($seqname) {
                if ($protWhole) {
                    $counterW++;
                    print_prot( "$tmpDir/$prot_addstop_file", $seqname,
                        $sequencepart );
                } else {
                    if ( defined( $protIDs{$seqname} ) ) {
                        for ( my $i = 0; $i < scalar( @{ $protIDs{$seqname} } ); $i++ ) {
                            if ($reg == 1) {
                                if ( defined( $contigIDs{ @{ $protIDs{$seqname} } [$i] } ) ) {
                                    $counterW++;
                                    print_prot("$tmpDir/@{$protIDs{$seqname}}[$i].fa", $seqname, $sequencepart);
                                }
                            } else {
                                $counterW++;
                                print_prot("$tmpDir/prot@{$protIDs{$seqname}}[$i].fa", $seqname, $sequencepart);
                            }
                        }
                    }
                }
            }
            $seqname = substr( $_, 1 );
            $sequencepart = "";
        } else {
            $sequencepart .= $_;
        }
    }

    if ($protWhole) {
        $counterW++;
        print_prot( "$tmpDir/$prot_addstop_file", $seqname, $sequencepart );
    } else {
        if ( defined( $protIDs{$seqname} ) ) {
            for ( my $i = 0; $i < scalar( @{ $protIDs{$seqname} } ); $i++ ) {
                if ($reg == 1) {
                    if (defined( $contigIDs{ @{ $protIDs{$seqname} }[$i] } ) ) {
                        print_prot( "$tmpDir/@{$protIDs{$seqname}}[$i].fa", $seqname, $sequencepart );
                    }
                } else {
                    print_prot( "$tmpDir/prot@{$protIDs{$seqname}}[$i].fa", $seqname, $sequencepart );
                }
            }
        }
    }
    close(PROT) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close file $prot_file!\n");
}

sub start_align {
    if ( !-d $alignDir ) {
        make_path($alignDir);
        print LOG "\# "
            . (localtime)
            . ": create working directory $alignDir\n";
        print LOG "mkdir $alignDir\n\n";
    } else {
        print LOG "\# "
            . (localtime)
            . ": working directory $alignDir already exists. Deleting files in that directory.\n";
        $cmdString = "rm -r $alignDir";
        system($cmdString) == 0 or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nFailed to execute $cmdString\n");
        make_path($alignDir);
    }

    chdir $tmpDir or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not change to directory $tmpDir.\n");

    my $pm                    = new Parallel::ForkManager($CPU);
    my $whole_prediction_file = "$alignDir/$prgsrc.concat.aln";
    if ($reg == 1) {
        foreach my $ID ( keys %contigIDs ) {
            my $pid = $pm->start and next;
            my $target = "$tmpDir/$contigIDs{$ID}{\"seq\"}$ID.fa";  # genome sequence
            my $query = "$ID.fa";                            # protein file
            $errorfile = "$alignDir/$contigIDs{$ID}{\"seq\"}.$ID.$prgsrc.stderr";
            $stdoutfile = "$alignDir/$contigIDs{$ID}{\"seq\"}.$ID.$prgsrc.aln";
            my $stdAdjusted;

            # create target file (contigIDs)
            if ( !-e $target ) {
                my $length = $contigIDs{$ID}{"end"} - $contigIDs{$ID}{"start"} + 1 + ( 2 * $offset );
                my $substart = $contigIDs{$ID}{"start"} - $offset;
                my $subseq   = substr( $seq{ $contigIDs{$ID}{"seq"} }, $substart, $length );
                print_seq( $target, $subseq, $contigIDs{$ID}{"seq"} );
            }
            if ( $prgsrc eq "exonerate" ) {
                call_exonerate( $target, $query, $stdoutfile, $errorfile );
                $stdAdjusted = adjust_exonerate( $stdoutfile, $ID );
            }

            if ( $prgsrc eq "spaln" ) {
                call_spaln( $target, $query, $stdoutfile, $errorfile, 1 );
                $stdAdjusted = adjust( $stdoutfile, $ID );
            }

            if ( $prgsrc eq "gth" ) {
                call_gth( $target, $query, $stdoutfile, $errorfile );
                $stdAdjusted = adjust( $stdoutfile, $ID );
            }
            $cmdString = "";
            if ($nice) {
                $cmdString .= "nice ";
            }
            $cmdString .= "cat $stdAdjusted >>$whole_prediction_file";
            print LOG "\# "
            . (localtime)
            . ": add prediction from file $stdAdjusted to file $whole_prediction_file\n";
            print LOG "$cmdString\n\n";
            system("$cmdString") == 0
            or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nfailed to execute: $cmdString!\n");
            unlink ( $stdAdjusted ) or die ("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nFailed to delete file $stdAdjusted!\n");
            $pm->finish;
        }
    } else {
        if ( $CPU > 1 || ( $CPU == 1 && !$protWhole ) ) {
            foreach my $ID ( keys %seq ) {
                my $pid = $pm->start and next;
                my $target = "$tmpDir/genome$ID.fa";    # genome sequence
                my $query;
                if ( !$protWhole ) {
                    $query = "prot$ID.fa";    # protein file
                }
                else {
                    $query = "$prot_addstop_file";
                }
                $errorfile  = "$alignDir/$ID.$prgsrc.stderr";
                $stdoutfile = "$alignDir/$ID.$prgsrc.aln";

                # create target file (whole sequences)
                if ( !-e $target ) {
                    print_seq( $target, $seq{$ID}, $ID );
                }
                if ( $prgsrc eq "exonerate" ) {
                    call_exonerate( $target, $query, $stdoutfile,
                        $errorfile );
                }

                if ( $prgsrc eq "spaln" ) {
                    if ( scalar( keys %seq ) > 1 ) {
                        call_spaln( $target, $query, $stdoutfile, $errorfile,
                            1 );
                    } else {
                        call_spaln( $target, $query, $stdoutfile, $errorfile,
                            $CPU );
                    }
                }

                if ( $prgsrc eq "gth" ) {
                    call_gth( $target, $query, $stdoutfile, $errorfile );
                }
                $cmdString = "";
                if ($nice) {
                    $cmdString .= "nice ";
                }
                $cmdString .= "cat $stdoutfile >>$whole_prediction_file";
                print LOG "\# "
                    . (localtime)
                    . ": add prediction from file $stdoutfile to file $whole_prediction_file\n";
                print LOG "$cmdString\n\n";
                system("$cmdString") == 0
                or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nfailed to execute: $cmdString!\n");
                unlink( $stdoutfile ) or die ("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nFailed to delete file $stdoutfile!\n");
                $pm->finish;
            }
        } else {
            $errorfile  = "$alignDir/$prgsrc.stderr";
            $stdoutfile = "$alignDir/$prgsrc.aln";
            if ( $CPU == 1 && $prgsrc eq "gth" && $protWhole ) {
                call_gth( $genome_file, "$prot_addstop_file", $stdoutfile,  $errorfile );
            } elsif ( $CPU == 1 && $prgsrc eq "exonerate" && $protWhole ) {
                call_exonerate( $genome_file, "$prot_addstop_file", $stdoutfile,  $errorfile);
            } elsif ( $CPU == 1 && $prgsrc eq "spaln" && $protWhole ) {
                call_spaln( $genome_file, "$tmpDir/$prot_addstop_file", $stdoutfile, $errorfile, 1 );
            }
            $cmdString = "";
            if ($nice) {
                $cmdString .= "nice ";
            }
            $cmdString .= "cat $stdoutfile >>$whole_prediction_file";
            print LOG "\# "
                . (localtime)
                . ": add prediction from file $stdoutfile to file $whole_prediction_file\n";
            print LOG "$cmdString\n\n";
            system("$cmdString") == 0
                or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nfailed to execute: $cmdString!\n");
            unlink( $stdoutfile ) or die ("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nFailed to delete $stdoutfile!\n");
        }
    }
    $pm->wait_all_children;
    chdir $dir or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not change to directory $dir.\n");

    #for version spaln2.2.0 and versions with same error
    if ( $prgsrc eq "spaln" ) {
        adjust_spaln_noreg($whole_prediction_file);
    }
}

# call assigned alignment program
sub call_exonerate {
    my $tFile      = shift;
    my $qFile      = shift;
    my $stdoutfile = shift;
    my $errorfile  = shift;
    $cmdString = "";
    if ($nice) {
        $cmdString .= "nice ";
    }
    if ( defined($ALIGNMENT_TOOL_PATH) ) {
        $cmdString .= $ALIGNMENT_TOOL_PATH . "/";
    }
    $cmdString .= "exonerate --model protein2genome --maxintron $maxintronlen --showtargetgff --showalignment no --query $qFile -t   $tFile ";
    if (defined($alnArgs)){
        $cmdString .= "$alnArgs ";
    }
    $cmdString .= "> $stdoutfile 2> $errorfile";
    print LOG "\# "
        . (localtime)
        . ": run exonerate for target '$tFile' and query '$qFile'\n";
    print LOG "$cmdString\n\n";
    system("$cmdString") == 0 or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nfailed to execute: $cmdString!\n");
}

sub call_gth {
    my $genomic = shift;
    my $protein = shift;

    # gth manipulates the protein file; if run in parallel, this leads to problems because multiple processes try to manipulate the same protein file
    # therefore we create a folder with a link to "protein", so that the linked file can be used by gth as template for file manipulations
    my @proteinPath = split( /\//, $protein );
    my $proteinStemPath;
    if ( scalar(@proteinPath) > 1 ) {
        for ( my $i = 0; $i < scalar(@proteinPath); $i++ ) {
            $proteinStemPath .= "/$proteinPath[$i]";
        }
    }
    my @genomePath = split( /\//, $genomic );
    $proteinStemPath .= "protFileFor_" . $genomePath[ ( scalar(@genomePath) - 1 ) ];
    mkdir $proteinStemPath;
    die "ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCannot create directory $proteinStemPath : $!\n" unless -d $proteinStemPath;
    $cmdString = "ln -s ../$protein $proteinStemPath/"
        . $proteinPath[ ( scalar(@proteinPath) - 1 ) ];
    print LOG "\# " . (localtime) . ": $cmdString\n";
    system("$cmdString") == 0 or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nfailed to execute: $cmdString!\n");
    $protein = $proteinStemPath . "/" . $proteinPath[ ( scalar(@proteinPath) - 1 ) ];
    my $stdoutfile = shift;
    my $errorfile  = shift;
    $cmdString = "";

    if ($nice) {
        $cmdString .= "nice ";
    }
    if ( defined($ALIGNMENT_TOOL_PATH) ) {
        $cmdString .= $ALIGNMENT_TOOL_PATH . "/";
    }
    if(not(defined($alnArgs))){
        $alnArgs = "-prseedlength $prseedlength -prminmatchlen $prminmatchlen -prhdist $prhdist";
    }
    $cmdString .= "gth -genomic $genomic -protein $tmpDir/$protein -gff3out -skipalignmentout -paralogs -gcmincoverage $gcmincoverage "
               . $alnArgs
               . " -o  $stdoutfile 2>$errorfile";
    print "CmdString: $cmdString\n";
    print LOG "\# " . (localtime) . ": run GenomeThreader for genome '$genomic' and protein '$protein'\n";
    print LOG "$cmdString\n\n";
    system("$cmdString") == 0 or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nfailed to execute: $cmdString!\n");
    $cmdString = "rm -r $proteinStemPath";
    print LOG "\# " . (localtime) . ": $cmdString\n";
    system("$cmdString") == 0 or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nfailed to execute: $cmdString!\n");
}

sub call_spaln {
    my $tFile      = shift;
    my $qFile      = shift;
    my $stdoutfile = shift;
    my $errorfile  = shift;
    my $cpus       = shift; # separate variable because multithreading shall only be used if data parallelization is not taking place, already
    my @split = split( /\./, $tFile );

    # some fo the spaln perl scripts expect spaln in $PATH
    $ENV{PATH} = "$ALIGNMENT_TOOL_PATH:$ENV{PATH}";
    if ( !-f "$split[0].idx" ) {
        $cmdString = "";
        if ($nice) {
            $cmdString .= "nice ";
        }
        if ( defined($ALIGNMENT_TOOL_PATH) ) {
            $cmdString .= $ALIGNMENT_TOOL_PATH . "/";
        }
        $cmdString .= "makdbs -KP $tFile ";
        print LOG "\# "
            . (localtime)
            . ": create *.ent, *.grp, *.idx, (*.odr), and *.seq files for target '$tFile' \n";
        print LOG "$cmdString\n\n";
        system("$cmdString") == 0 or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nfailed to execute: $cmdString!\n");
    }
    if ( !-f "$split[0].bkp" ) {
        $cmdString = "";
        if ($nice) {
            $cmdString .= "nice ";
        }
        $cmdString .= "perl $ENV{'ALN_DBS'}/makblk.pl -W$split[0].bkp -KP $tFile";
        print LOG "\# "
            . (localtime)
            . ": create *.bkp file for target '$tFile'\n";
        print LOG "$cmdString\n\n";
        system("$cmdString") == 0 or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nfailed to execute: $cmdString!\n");
    }
    @split = split( /\./, $qFile );
    if ( !-f "$split[0].idx" ) {
        $cmdString = "";
        $cmdString = "";
        if ($nice) {
            $cmdString .= "nice ";
        }
        if ( defined($ALIGNMENT_TOOL_PATH) ) {
            $cmdString .= $ALIGNMENT_TOOL_PATH . "/";
        }
        $cmdString .= "makdbs -KA $qFile";
        print LOG "\# " . (localtime)
            . ": create *.ent, *.grp, *.idx, (*.odr), and *.seq files for query '$qFile' \n";
        print LOG "$cmdString\n\n";
        system("$cmdString") == 0 or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nfailed to execute: $cmdString!\n");
    }
    if ( !-f "$split[0].bka" ) {
        $cmdString = "";
        if ($nice) {
            $cmdString .= "nice ";
        }
        $cmdString .= "perl $ENV{'ALN_DBS'}/makblk.pl -W$split[0].bka -KA $qFile";
        print LOG "\# " . (localtime)
            . ": create *.bka file for query '$qFile'\n";
        print LOG "$cmdString\n\n";
        print LOG "\# OUTPUT OF makblk.pl STARTING #\n";
        system("$cmdString") == 0 or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nfailed to execute: $cmdString!\n");
    }
    $cmdString = "";
    if ($nice) {
        $cmdString .= "nice ";
    }
    if ( defined($ALIGNMENT_TOOL_PATH) ) {
        $cmdString .= $ALIGNMENT_TOOL_PATH . "/";
    }
    $cmdString .= "spaln ";
    if ( $cpus > 1 ) {
        $cmdString .= "-t[$cpus] ";
    }
    $cmdString .= "-Q7 -O0 ";
    if(defined($alnArgs)){
        $cmdString .=  "$alnArgs ";
    }
    $cmdString .= "$tFile $qFile > $stdoutfile 2>$errorfile";
    print LOG "\# " . (localtime)
        . ": run spaln for target '$tFile' and query '$qFile'\n";
    print LOG "$cmdString\n\n";
    system("$cmdString") == 0 or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nfailed to execute: $cmdString!\n");
}

sub print_prot {
    my $outfile = shift;
    my $Sname   = shift;
    my $seqpart = shift;
    open( OUT, ">>$outfile" ) or die "ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCannot open file: $outfile\n";
    print OUT ">$Sname\n";

    # add '*' for GenomeThreader (gth) [instead of 'gt seqtransform -addstopaminos ...' since this is not part
    # of gth, but of GenomeTools (gt)]
    if ( $prgsrc eq "gth" ) {
        print OUT "$seqpart*\n";
        if ( $counterW == 1 ) {
            print LOG "\# " . (localtime)
                . ": add stop amino acid ('*') to protein sequences'\n";
        }
    } else {
        print OUT "$seqpart\n";
    }
    close(OUT) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close file $outfile!\n");
}

# print genome sequence part or whole sequence
sub print_seq {
    my $tFile    = shift;
    my $sequence = shift;
    my $Sname    = shift;
    open( TARGET, ">$tFile" ) or die "ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCannot open file: $tFile\n";
    print TARGET ">$Sname\n";
    my $start = 0;     # see getAnnoFasta.pl
    my $ret   = "";    # see getAnnoFasta.pl
    while ( length($sequence) - $start >= 80 ) {    # see getAnnoFasta.pl
        my $shortline = substr( $sequence, $start, 80 );
        $ret .= "$shortline\n";
        $start += 80;
    }
    $ret .= substr( $sequence, $start, 80 ) . "\n"
    if ( $start < length($sequence) );
    print TARGET $ret;
    close(TARGET) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close file $tFile!\n");
}

# adjust spaln and gth gff3 output
sub adjust {
    my $output     = shift;
    my $conID      = shift;
    my $output_adj = "$output.adj";
    open( OUT, ">$output_adj" ) or die "ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCannot open file: $output_adj\n";
    open( IN,  $output )        or die "ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCannot open file: $output\n";
    while (<IN>) {
        chomp;
        my @f = split( /\s/, $_ );
        if ( $_ =~ m/##sequence-region/ ) {
            $f[-1] += $contigIDs{$conID}{"start"} - $offset + $spalnErrAdj;
            $f[-2] += $contigIDs{$conID}{"start"} - $offset + $spalnErrAdj;
            print OUT "##sequence-region\t$f[-3] $f[-2] $f[-1]\n";
        } elsif ( scalar(@f) >= 9 ) {
            $f[3] += $contigIDs{$conID}{"start"} - $offset + $spalnErrAdj;
            $f[4] += $contigIDs{$conID}{"start"} - $offset + $spalnErrAdj;
            print OUT
            "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\t$f[8]\n";
        } else {
            print OUT "$_\n";
        }
    }
    close(IN)  or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close file $output!\n");
    close(OUT) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close file $output_adj!\n");
    return $output_adj;
}

# adjust exonerate output
sub adjust_exonerate {
    my $exonerate_output     = shift;
    my $conID                = shift;
    my $exonerate_output_adj = "$exonerate_output.adj";
    open( OUT, ">$exonerate_output_adj" )
    or die "ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCannot open file: $exonerate_output_adj\n";
    open( IN, $exonerate_output )
    or die "ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCannot open file: $exonerate_output\n";
    my $add_space = "";
    while (<IN>) {
        chomp;
        if ( $_ =~ m/  Target range: (\d+) -> (\d+)/ ) {
            my $start = $1 + $contigIDs{$conID}{"start"} - $offset;
            my $end   = $2 + $contigIDs{$conID}{"start"} - $offset;
            my $diff  = length($end) - length($2);
            $add_space = " " x $diff;
            print OUT "  Target range: $start -> $end\n";
        } elsif ( $_ =~ m/^(\s+)(\d+) : ([^a-z\.\-]+) : (\d+)/ ) {
            my $start = $2 + $contigIDs{$conID}{"start"} - $offset;
            my $end   = $4 + $contigIDs{$conID}{"start"} - $offset;
            print OUT "$1" . "$start : $3 : " . "$end\n";
        } elsif ( $_ =~ m/vulgar:/ ) {
            my @line = split( /[\s, ]/, $_ );
            $line[6] += $contigIDs{$conID}{"start"} - $offset;
            $line[7] += $contigIDs{$conID}{"start"} - $offset;
            my $vulgar = join( ' ', @line );
            print OUT "$vulgar\n";
        } elsif ($_ =~ m/^(\s+)(\d+) : ([a-zA-Z\.\-]+)/
            || $_ =~ m/^(\s+)([a-zA-Z\.\-]+)/
            || $_ =~ m/^(\s+)[\|\.\!\:]+/ )
        {
            print OUT $add_space . $_ . "\n";
        } else {
            my @f = split( /\t/, $_ );
            if ( scalar(@f) == 9 || scalar(@f) == 8 ) {
                if ( $f[6] eq "+" ) {
                    $f[3] += $contigIDs{$conID}{"start"} - $offset;
                    $f[4] += $contigIDs{$conID}{"start"} - $offset;
                }
                else {
                    $f[3] += $contigIDs{$conID}{"start"} - $offset;
                    $f[4] += $contigIDs{$conID}{"start"} - $offset;
                }
                print OUT
                "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]";
                if ( scalar(@f) == 9 ) {
                    if ( $f[2] eq "similarity" ) {
                        my @sim_f = split( /\t/, $f[8] );
                        for ( my $i = 0; $i < scalar(@sim_f); $i++ ) {
                            if ( $sim_f[$i] eq "Align" ) {
                                $sim_f[ $i + 1 ]
                                += $contigIDs{$conID}{"start"};
                            }
                        }
                        $f[8] = join( ' ', @sim_f );
                    }
                    print OUT "\t$f[8]";
                }
                print OUT "\n";
            } else {
                print OUT "$_\n";
            }
        }
    }
    close(IN)  or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close file $exonerate_output!\n");
    close(OUT) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close file $exonerate_output_adj!\n");
    return $exonerate_output_adj;
}

# adjust spaln output for noreg option
sub adjust_spaln_noreg {
    my $output     = shift;
    my $output_adj = "$output.adj";
    open( OUT, ">$output_adj" ) or die "ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCannot open file: $output_adj\n";
    open( IN,  $output )        or die "ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCannot open file: $output\n";
    while (<IN>) {
        chomp;
        my @f = split( /\s/, $_ );
        if ( $_ =~ m/##sequence-region/ ) {
            $f[-1] += $spalnErrAdj;
            $f[-2] += $spalnErrAdj;
            print OUT "##sequence-region\t$f[-3] $f[-2] $f[-1]\n";
        } elsif ( scalar(@f) >= 9 ) {
            $f[3] += $spalnErrAdj;
            $f[4] += $spalnErrAdj;
            print OUT
            "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\t$f[8]\n";
        } else {
            print OUT "$_\n";
        }
    }
    close(IN)  or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close file $output!\n");
    close(OUT) or die("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nCould not close file $output_adj!\n");
}

# delete empty files
sub clean_up {
    print LOG "NEXT STEP: delete empty files\n";
    my @files = `find $alignDir -empty`;
    print LOG "\# " . (localtime) . ": delete empty files\n";
    for ( my $i = 0; $i <= $#files; $i++ ) {
        # to prevent error: Unsuccessful stat on filename containing newline
        chomp( $files[$i] );
        if ( -f $files[$i] ) {
            print LOG "rm $files[$i]\n";
            unlink( rel2abs( $files[$i] ) );
        }
    }
    print LOG "empty files deleted.\n";
    print LOG "NEXT STEP: deleting $tmpDir\n";
    rmtree( ["$tmpDir"] ) or die ("ERROR in file " . __FILE__ ." at line ". __LINE__ ."\nFailed to delete folder $tmpDir!\n");
    print LOG "directory deleted.\n";
}
