####################################################################################################
#                                                                                                  #
# helpMod_braker.pm                                                                                       #
#                                                                                                  #
# Component of braker.pl                                                                           #
#                                                                                                  #
# Authors: Katharina Hoff, Simone Lange, Mario Stanke, Alexandre Lomsadze, Mark Borodovsky         #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #                                                             #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
#                                                                                                  #
####################################################################################################

package helpMod_braker;

use Exporter 'import';

@EXPORT_OK = qw( find tildeConvert checkFile formatDetector relToAbs setParInConfig addParToConfig
    uptodate gtf2fasta clean_abort );

use strict;
use Cwd;
use Cwd 'abs_path';
use File::Spec::Functions qw(rel2abs);
use File::Basename qw(dirname);
use File::Path qw(rmtree);

####################################################################################################
# genetic code (use only one for training gene blast all-against-all, not gene prediction)         #
####################################################################################################
my %genetic_code;
$genetic_code{1} = {
    'TCA' => 'S', # Serine
    'TCC' => 'S', # Serine
    'TCG' => 'S', # Serine
    'TCT' => 'S', # Serine
    'TTC' => 'F', # Phenylalanine
    'TTT' => 'F', # Phenylalanine
    'TTA' => 'L', # Leucine
    'TTG' => 'L', # Leucine
    'TAC' => 'Y', # Tyrosine
    'TAT' => 'Y', # Tyrosine
    'TAA' => '*', # Stop
    'TAG' => '*', # Stop
    'TGC' => 'C', # Cysteine
    'TGT' => 'C', # Cysteine
    'TGA' => '*', # Stop
    'TGG' => 'W', # Tryptophan
    'CTA' => 'L', # Leucine
    'CTC' => 'L', # Leucine
    'CTG' => 'L', # Leucine
    'CTT' => 'L', # Leucine
    'CCA' => 'P', # Proline
    'CAT' => 'H', # Histidine
    'CAA' => 'Q', # Glutamine
    'CAG' => 'Q', # Glutamine
    'CGA' => 'R', # Arginine
    'CGC' => 'R', # Arginine
    'CGG' => 'R', # Arginine
    'CGT' => 'R', # Arginine
    'ATA' => 'I', # Isoleucine
    'ATC' => 'I', # Isoleucine
    'ATT' => 'I', # Isoleucine
    'ATG' => 'M', # Methionine
    'ACA' => 'T', # Threonine
    'ACC' => 'T', # Threonine
    'ACG' => 'T', # Threonine
    'ACT' => 'T', # Threonine
    'AAC' => 'N', # Asparagine
    'AAT' => 'N', # Asparagine
    'AAA' => 'K', # Lysine
    'AAG' => 'K', # Lysine
    'AGC' => 'S', # Serine
    'AGT' => 'S', # Serine
    'AGA' => 'R', # Arginine
    'AGG' => 'R', # Arginine
    'CCC' => 'P', # Proline
    'CCG' => 'P', # Proline
    'CCT' => 'P', # Proline
    'CAC' => 'H', # Histidine
    'GTA' => 'V', # Valine
    'GTC' => 'V', # Valine
    'GTG' => 'V', # Valine
    'GTT' => 'V', # Valine
    'GCA' => 'A', # Alanine
    'GCC' => 'A', # Alanine
    'GCG' => 'A', # Alanine
    'GCT' => 'A', # Alanine
    'GAC' => 'D', # Aspartic Acid
    'GAT' => 'D', # Aspartic Acid
    'GAA' => 'E', # Glutamic Acid
    'GAG' => 'E', # Glutamic Acid
    'GGA' => 'G', # Glycine
    'GGC' => 'G', # Glycine
    'GGG' => 'G', # Glycine
    'GGT' => 'G'  # Glycine
};

$genetic_code{6} = { %{$genetic_code{1}} };
$genetic_code{6}{'TAA'} = 'Q'; # Glutamine
$genetic_code{6}{'TAG'} = 'Q';  # Glutamine

$genetic_code{10} = { %{$genetic_code{1}} };
$genetic_code{10}{'TGA'} = 'C'; # Cysteine

$genetic_code{12} = { %{$genetic_code{1}} };
$genetic_code{12}{'CTG'} = 'S'; # Serine

$genetic_code{25} = { %{$genetic_code{1}} };
$genetic_code{25}{'TGA'} = 'G'; # Glycine

$genetic_code{26} = { %{$genetic_code{1}} };
$genetic_code{26}{'CTG'} = 'A'; # Alanine

$genetic_code{27} = { %{$genetic_code{1}} };
$genetic_code{27}{'TAG'} = 'Q'; # Glutamine
$genetic_code{27}{'TAA'} = 'Q';  # Glutamine
# cannot differentiate between TGA Stop or Tryptophane, make stop codon always

$genetic_code{28} = { %{$genetic_code{1}} };
# cannot differentiate between alternative translation of stop codons, keep table 1

$genetic_code{29} = { %{$genetic_code{1}} };
$genetic_code{29}{'TAA'} = 'Y'; # Tyrosine
$genetic_code{29}{'TAG'} = 'Y';  # Tyrosine

$genetic_code{30} = { %{$genetic_code{1}} };
$genetic_code{30}{'TAA'} = 'E'; # Glutamic Acid
$genetic_code{30}{'TAG'} = 'E';  # Glutamic Acid

$genetic_code{31} = { %{$genetic_code{1}} };
$genetic_code{31}{'TGA'} = 'W'; # Tryptophane
# cannot differentiale alternative translation of other codons, keep table 1



####################################################################################################
# extract DNA sequence of CDS in gtf from genome fasta file, write to CDS fasta file               #
####################################################################################################
sub gtf2fasta {
    my $genome_file = shift;
    my $gtf_file = shift;
    my $fasta_file = shift;
    my $table = shift;
    my %gtf;
    my %genome;
    open (GTF, "<", $gtf_file ) or die ("ERROR: in file " . __FILE__ ." at line ". __LINE__
        . "\nCould not close file $gtf_file!\n");
    while ( <GTF> ) {
        if ( $_ =~ m/\tCDS\t/ ) {
            $_ =~ m/transcript_id \"(\S+)\"/;
            my $txid = $1;
            my @line = split(/\t/);
            push @{$gtf{$line[0]}{$txid}}, $_;
        }
    }
    close (GTF) or die ("ERROR: in file " . __FILE__ ." at line ". __LINE__
        . "\nCould not close file $gtf_file!\n");
    open (GENOME, "<", $genome_file ) or die ("ERROR: in file " . __FILE__
        . " at line ". __LINE__ ."\nCould not close file $genome_file!\n");
    my $seq = "";
    my $locus;
    my %cds_seq;
    while (<GENOME>) {
        chomp;
        if( not ($_ =~ m/^>/ ) ) {
            $seq .= $_;
        } elsif ( $_ =~ m/^>(\S+)/ ) {
            if (defined ($locus) ) {
                while ( my ( $txid, $txgtf ) = each %{$gtf{$locus}} ) {
                    foreach ( @{$txgtf} ) {
                        my @line = split(/\t/);
                        if ( not( defined( $cds_seq{$txid} ) ) ) {
                            $cds_seq{$txid} = substr ( $seq, ( $line[3] -1 ),
                                ( $line[4] - $line[3] + 1 ) );
                            if ( $line[6] eq '-' ) {
                                $cds_seq{$txid} = reverse_complement ( $cds_seq{$txid} );
                            }
                        }else {
                            if ( $line[6] eq '+') {
                                $cds_seq{$txid} .= substr ( $seq, ( $line[3] -1 ),
                                    ( $line[4] - $line[3] + 1 ) );
                            } else {
                                $cds_seq{$txid} = reverse_complement(
                                    substr ( $seq, ( $line[3] -1 ), ( $line[4] - $line[3] + 1 ) ) )
                                    . $cds_seq{$txid};
                            }
                        }
                    }
                }
            }
            $locus = $1;
            $seq = "";
        }
    }
    # excise seqs for last contig:
    if (defined ($locus) ) {
        while ( my ( $txid, $txgtf ) = each %{$gtf{$locus}} ) {
            foreach ( @{$txgtf} ) {
                my @line = split(/\t/);
                if ( not( defined( $cds_seq{$txid} ) ) ) {
                    $cds_seq{$txid} = substr ( $seq, ( $line[3] -1 ), ( $line[4] - $line[3] + 1 ) );
                    if ( $line[6] eq '-' ) {
                        $cds_seq{$txid} = reverse_complement ( $cds_seq{$txid} );
                    }
                }else {
                    if ( $line[6] eq '+') {
                        $cds_seq{$txid} .= substr ( $seq, ( $line[3] -1 ),
                            ( $line[4] - $line[3] + 1 ) );
                    } else {
                        $cds_seq{$txid} = reverse_complement(
                            substr ( $seq, ( $line[3] -1 ), ( $line[4] - $line[3] + 1 ) ) )
                            . $cds_seq{$txid}
                    }
                }
            }
        }
    }
    close(GENOME) or die ("ERROR: in file " . __FILE__ ." at line ". __LINE__
        . "\nCould not close file $genome_file!\n");
    open (FASTA, ">", $fasta_file) or die ("ERROR: in file " . __FILE__ ." at line " . __LINE__
        . "\nCould not close file $fasta_file!\n");
    while ( my ( $txid, $dna ) = each %cds_seq ) {
        print FASTA ">$txid\n".dna2aa($dna, $table)."\n";
    }
    close (FASTA) or die ("ERROR: in file " . __FILE__ ." at line ". __LINE__
        . "\nCould not close file $fasta_file!\n");
}

####################################################################################################
# reverse complement of DNA sequence                                                               #
####################################################################################################

sub reverse_complement {
    my $in = shift;
    $in =~ tr/ACGTacgt/TGCAtgca/;
    $in = reverse ($in);
    return $in;
}

####################################################################################################
# Translate DNA to protein sequence                                                                #
####################################################################################################

sub dna2aa {
    my $seq = shift;
    my $code = shift;
    $seq = uc($seq);
    my @codons = $seq =~ /(.{1,3})/g;
    my $aa = "";
    foreach ( @codons ) {
        if($_ =~ m/N/i){
            $aa .= "X";
        }else{
            $aa .= $genetic_code{$code}{$_};
        }
    }
    return $aa;
}

####################################################################################################
# find a script a script that is used by braker.pl                                                 #
####################################################################################################

sub find {
    my $script                = shift;    # script to find
    my $augustus_bin_path     = shift;    # augustus_bin_path
    my $augustus_scripts_path = shift;    # augustus_scripts_path
    my $AUGUSTUS_CONFIG_PATH  = shift;
    my $exist;     # boolean variable, to remark if $script is found
    my $string;    # the absolute name for $script

    my $path_1
        = abs_path($augustus_scripts_path);  # first searching path of $script
    my $path_2 = abs_path("$augustus_bin_path/../scripts")
        ;    # second searching path of $script
    my $path_3 = abs_path("$AUGUSTUS_CONFIG_PATH/../scripts")
        ;    # third searching path of $script
    my $path_4 = dirname( rel2abs(__FILE__) );

    # paths can be redundant, remove redundancies:
    my @unique;
    my %seen;
    foreach my $value ( ( "$path_1/$script", "$path_2/$script", "$path_3/$script",
        "$path_4/$script" ) )
    {
        if ( !$seen{$value}++ ) {
            push @unique, $value;
        }
    }
    foreach (@unique) {
        if ( -f $_ ) {
            $exist  = 1;
            $string = $_;
            last;
        }
    }
    if ($exist) {

        # print "Found script $string.\n";
        return $string;    # return name with absolute path
    }
    else {
        # if not found, output error
        die("ERROR: in file " . __FILE__ ." at line ". __LINE__
            . "\n found neither $path_1/$script nor $path_2/$script nor $path_3/$script nor "
            . "$path_4/$script!\nPlease Check the environment variables AUGUSTUS_CONFIG_PATH and "
            . "command line options AUGUSTUS_BIN_PATH and AUGUSTUS_SCRIPTS_PATH or install "
            . "AUGUSTUS, again!\n"
        );
    }
}

####################################################################################################
# convert a file name which begins with ~ to name with absolute path                               #
####################################################################################################

sub tildeConvert {
    my $file = shift;
    if ( $file =~ /^~/ ) {
        my $HOME = $ENV{'HOME'};
        $file =~ s/~/$HOME/;    # replace ~ with home directory
    }
    return $file;
}

####################################################################################################
# check if $file exists and replace $file with its absolute path                                   #
####################################################################################################

sub checkFile {
    my $file = shift;           # file which to be checked
    my $type = shift
        ;   # type of file, used by error outputting if the file doesn't exist
    my $usage = shift;    # usage to be outputted if the file doesn't exist

    die("ERROR: in file " . __FILE__ ." at line ". __LINE__
        . "\nmissing $type file!\n$usage") if ( !$file );

    # overwrite $file with absolute path
    $file = tildeConvert($file);
    $file = rel2abs($file);        # overwrite $file with absolute path
    if ( !( -f $file ) ) {
        die("ERROR: in file " . __FILE__ ." at line ". __LINE__
            . "\n $type file $file not found!\n");
    }
    return $file;
}

####################################################################################################
# detect if a file has gff or gb or fasta format                                                   #
####################################################################################################

sub formatDetector {
    my $file      = shift;         # file to be detected
    my $testLines = 1000;          # read at most this many lines for testing
    my $i;
    my @helpArray_gff;
    #
    # check if file has GENBANK format
    #
    open( DFILE, $file ) or die("ERROR: in file " . __FILE__ ." at line ". __LINE__ .
        "\nCould not open $file!\n");
    $i = 0;
    my $haveLOCUS    = 0;
    my $haveSource   = 0;
    my $haveOrigin   = 0;
    my $haveTermSymb = 0;
    while ( defined( my $line = <DFILE> ) && $i < $testLines ) {
        $i++;
        $haveLOCUS++ if ( $i == 1 && $line =~ m/^LOCUS/ );
        $haveSource++   if ( $line =~ m/ +source +/i );
        $haveOrigin++   if ( $line =~ m/^ORIGIN/ );
        $haveTermSymb++ if ( $line =~ m/\/\// );
    }
    close(DFILE);
    if ((     ( $haveLOCUS > 0 )
            + ( $haveSource > 0 )
            + ( $haveOrigin > 0 )
            + ( $haveTermSymb > 0 )
        ) > 1
        )
    {
        print STDERR
            "ERROR: in file " . __FILE__ ." at line ". __LINE__
            . "\n$file appears to be in corrupt Genbank format. 'LOCUS' missing\n"
            if ( !$haveLOCUS );
        print STDERR
            "ERROR: in file " . __FILE__ ." at line ". __LINE__
            . "\n$file appears to be in corrupt Genbank format. ' source ' line missing\n"
            if ( !$haveSource );
        print STDERR
            "$file appears to be in corrupt Genbank format. 'ORIGIN' missing\n"
            if ( !$haveOrigin );
        print STDERR
            "ERROR: in file " . __FILE__ ." at line ". __LINE__
            . "\n$file appears to be in corrupt Genbank format. '//' missing\n"
            if ( !$haveTermSymb );
        return "gb";
    }
    #
    # check if file has GFF format
    #
    open( DFILE, $file ) or die("ERROR: in file " . __FILE__ ." at line ". __LINE__
        . "\nCould not open $file!\n");
    $i = 0;
    my $badGFFlines  = 0;
    my $goodGFFlines = 0;
    while ( defined( my $line = <DFILE> ) && $i < $testLines ) {
        $i++;

      # if not genbank format and the row not a possible comment in gff format
        if ( !( $line =~ m/^#/ ) && !( $line =~ m/^\s*$/ ) ) {
            @helpArray_gff = split( /\t/, $line );
            if ( $#helpArray_gff < 7 ) {

# each non-comment row should contain at least 7 tabulators (end with new line???)
                $badGFFlines++;
            }
            else {
                $goodGFFlines++;
            }
        }
    }
    close(DFILE);
    if ( $goodGFFlines > 0 ) {
        if ( $badGFFlines > 0 ) {
            print STDERR "ERROR: in file " . __FILE__ ." at line ". __LINE__
            . "\n$file appears to be in corrupt GFF format.\n";
            return "";
        }
        else {
            return "gff";
        }
    }
    #
    # check if file has FASTA format and whether it is DNA or protein
    #
    open( DFILE, $file ) or die("ERROR: in file " . __FILE__ ." at line ". __LINE__
        . "\nCould not open $file!\n");
    $i = 0;
    my $greaterLines = 0;
    my $concatseq    = "";
    while ( defined( my $line = <DFILE> ) && $i < $testLines ) {
        $i++;
        chomp $line;
        if ( $line =~ m/^>/ ) {
            $greaterLines++;
        }
        else {
            $concatseq .= $line;
        }
    }
    close(DFILE);
    if ( $greaterLines > 0 ) {
        $concatseq = uc($concatseq);
        my $len      = length($concatseq);
        my $protchar = ( $concatseq =~ tr/ACDEFGHIKLMNPQRSTVWYX// );
        my $dnachar  = ( $concatseq =~ tr/ACGTN// );
        if ( $protchar > 0.8 * $len ) {
            return "fasta-prot";
        }
        elsif ( $dnachar > 0.8 * $len ) {
            return "fasta-dna";
        }
        else {
            print STDERR "ERROR: in file " . __FILE__ ." at line ". __LINE__
            . "\n$file appears to be in corrupt FASTA format.\n";
            return "";
        }
    }
    return "";
}

####################################################################################################
# convert relative path to absolute path                                                           #
####################################################################################################

sub relToAbs {
    my $name = shift;
    $name = tildeConvert($name);    # overwrite working directory
    return rel2abs($name);          # with absolute path
}

####################################################################################################
# change a parameter in a config file                                                              #
# assume the format                                                                                #
# parName    value   # comment                                                                     #
####################################################################################################

sub setParInConfig {
    my $configFileName = shift;
    my $parName        = shift;
    my $value          = shift;
    open( CFGFILE, "+<$configFileName" )
        or die("Could not read config file $configFileName\n");
    my @lines = <CFGFILE>;
    foreach my $line (@lines) {
        $line =~ s/(\s*$parName\s+)(\S+?)(\s|\#|$)/$1$value$3/;
    }
    seek( CFGFILE, 0, 0 );
    print CFGFILE @lines or die("Could not write $configFileName");
    truncate( CFGFILE, tell(CFGFILE) );
    close(CFGFILE);
}

####################################################################################################
# add a parameter in to a config file                                                              #
# assume the format                                                                                #
# parName    value   # comment                                                                     #
####################################################################################################

sub addParToConfig {
    my $configFileName = shift;
    my $parName        = shift;
    my $value          = shift;
    open( CFGFILE, "+<$configFileName" )
        or die("Could not read config file $configFileName\n");
    my @lines = <CFGFILE>;
    push(@lines, '# ADDITIONAL PARAMTERS ADDED TO CFG BY BRAKER:\n');
    push(@lines, '$parName    $value\n');
    seek( CFGFILE, 0, 0 );
    print CFGFILE @lines or die("Could not write $configFileName");
    truncate( CFGFILE, tell(CFGFILE) );
    close(CFGFILE);
}

####################################################################################################
# uptodate                                                                                         #
# check whether output files are up to date with respect to input files                            #
# all output files must exist and not be older than any input file                                 #
####################################################################################################

sub uptodate {
    my $input  = shift;    # reference to list of input file names
    my $output = shift;    # reference to list of output file names
    my $earliestOutMtime;  # earliest modification time of an output file
    my $latestInMtime;     # latest modification time an any input file
    my @stat;              # holds info about file
    return 1 if ( @{$output} == 0 );    # no output is always up to date
                                        # check whether all output files exist
    foreach my $of ( @{$output} ) {
        return 0 if ( !-f $of );
        @stat             = stat($of);
        $earliestOutMtime = $stat[9]
            if ( !defined($earliestOutMtime)
            || $stat[9] < $earliestOutMtime );
    }
    return 1
        if ( @{$input} == 0 );    # no input is always older than output files
                                  # check existence and times of input files
    foreach my $if ( @{$input} ) {
        if ( !-f $if ) {          # ignore if input file does not exist
            print STDERR "Warning: $if missing.\n"
                ;                 # TODO, remove or correct this later
        }
        @stat          = stat($if);
        $latestInMtime = $stat[9]
            if ( !defined($latestInMtime) || $stat[9] > $latestInMtime );
    }
    return ( $latestInMtime <= $earliestOutMtime );
}

####################################################################################################
# exit braker after deleting AUGUSTUS parameter directory                                          #
# use instead of exit(1) and die() after creating of                                               #
# directory and before first real etraining                                                        #
####################################################################################################
sub clean_abort {
    my $configDir = shift;
    my $useexisting = shift;
    my $message = shift;
    if (-d $configDir && not($useexisting)) {
        rmtree( ["$configDir"] ) or die ("ERROR in file " . __FILE__ ." at line ". __LINE__
            . "\nFailed to delete $configDir!\n");
    }
    print STDERR $message;
    exit(1);
}

1;
