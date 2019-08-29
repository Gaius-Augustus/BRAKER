#!/usr/bin/env perl

####################################################################################################
#                                                                                                  #
# align2hints.pl - generate hints from spaln [O0 (=gff3)], exonerate, genomeThreader (gth)         #
#                  or scipio output                                                                #
#                                                                                                  #
# Authors: Katharina Hoff, Simone Lange, Mario Stanke                                              #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# Last modification: August 29th 2019                                                              #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
#                                                                                                  #
# Usage:                                                                                           #
# align2hints.pl [OPTIONS] --in=align.gff3 --out=hintsfile.gff --prg=gth|exonerate|spaln|scipio    #
#                                                                                                  #
####################################################################################################

use strict;
use warnings;
use Getopt::Long;
use File::Spec::Functions qw(rel2abs);
use Cwd;

my $usage = <<'ENDUSAGE';

align2hints.pl    generate hints from spaln [O0 (=gff3)], exonerate,
                  GenomeThreader (gth), scipio
                  or GEMOMA output.
                  Spaln2 run like this: spaln -O0 ... > spalnfile
                  Exonerate run like this:
                      exonerate --model protein2genome --showtargetgff T \
                         ... > exfile
                  GenomeThreader run like this: 
                      gth -genomic genome.fa  -protein protein.fa -gff3out \
                         -skipalignmentout ... -o gthfile
                  scipio run like this:
                  scipio.1.4.1.pl genome.fa prot.fa | yaml2gff.1.4.pl \
                      > scipio.gff

SYNOPSIS

align2hints.pl [OPTIONS] --in=align.gff3 --out=hintsfile.gff \
                         --prg=gth|exonerate|spaln|scipio

  --in                 input file from gth (gff3), spaln (gff3) or exonerate
                       output
  --out                contains CDSpart, CDS and intron hints


OPTIONS

    --help                   Print this help message.
    --CDSpart_cutoff=n       This many bp are cut off of each CDSpart hint
                             w.r.t. the cds (default 15).
    --maxintronlen=n         Alignments with longer gaps are discarded
                             (default 350000).
    --minintronlen=n         Alignments with gaps shorter than this and longer
                             than maxgaplen are discarded (default 41).
    --priority=n             Priority of hint group (default 4).
    --prg=s                  Alignment program of input file, either 'gth',
                             'spaln', 'exonerate', 'scipio', or 'gemoma'.
    --source=s               Source identifier (default 'P')
    --CDS                    Do not output CDSpart hints, but complete CDS
                             hints.
    --genome_file=s          if prg is exonerate and start hints shall be
                             created, the genome file from which the
                             alignments were generated, must be specified.
    --version                print version of align2hints.pl

Format:
  seqname <TAB> source <TAB> feature <TAB> start <TAB> end <TAB> score <TAB>
     strand <TAB> frame <TAB> src=source;grp=target_protein;pri=priority


DESCRIPTION

  Example:

    align2hints.pl [OPTIONS] --in=align.gff3 --out=hintsfile.gff --prg=gth

ENDUSAGE

my $version = 1.0;    # version of align2hints.pl
my $printVersion;
my $alignfile;        # alignment input file
my $CDSpart_cutoff = 15;           # cutoff for CDSpart hints
my $CDSpartid      = "CDSpart";    # abbreviate to decrease file size
my $dir;               # working superdirectory where programme is called from
my $hintsfilename;     # hints file output name
my $intron_end;        # end of intron (only for gth and spaln)
my $intron_score;      # intron score
my $intron_start;      # start of intron (only for gth and spaln)
my $intron_threshold;  # Threshold for current programme
my $intron_threshold_gth
    = 0.7;    # introns from gth with a score below this are discarded
my $intron_threshold_spn
    = 200;    # introns from spaln with a score below this are discarded
my $maxintronlen = 350000;    # maximal intron length
my $minintronlen = 41;        # default minimal intron length
my $parent;                   # current parent
my $prevParent = "noP"; # previous parent
my $prevScore  = 0;     # previous exon/CDS score for calculating intron score
# positions in query protein (scipio only), to determine intron for scipio
my ( $qstart, $qend, $prevQend );
my $prgsrc;    # source programme (exonerate, spaln or gth)
my $priority = 4;      # priority for hints
my $source   = "P";    # source for extrinsic file
my $genome_file;       # genome file name
my %genome;            # hash to store genome sequence
my $CDS; # output CDS instead of CDSpart hints
my $help;

if ( $#ARGV < 1 || $help ) {
    print $usage;
    exit;
}

GetOptions(
    'in=s'             => \$alignfile,
    'out=s'            => \$hintsfilename,
    'CDSpart_cutoff:i' => \$CDSpart_cutoff,
    'dir=s'            => \$dir,
    'maxintronlen:i'   => \$maxintronlen,
    'minintronlen:i'   => \$minintronlen,
    'priority:i'       => \$priority,
    'prg=s'            => \$prgsrc,
    'genome_file=s'    => \$genome_file,
    'CDS!'             => \$CDS,
    'help!'            => \$help,
    'source:s'         => \$source,
    'version!'         => \$printVersion
);

if ($printVersion) {
    print "align2hints.pl version $version\n";
    exit(0);
}

# mainly for usage within BRAKER2
if ( !defined($dir) ) {
    $dir = cwd();
}
my $last_char = substr( $dir, -1 );
if ( $last_char eq "\/" ) {
    chop($dir);
}

if ( defined($alignfile) ) {

    # check whether alignment file exists
    if ( !-e $alignfile ) {
        print STDERR
            "ERROR: Alignment file $alignfile does not exist. Please check.\n";
        exit(1);
    }
    else {
        $alignfile = rel2abs($alignfile);
    }
}

if ( !defined($prgsrc) ) {
    print STDERR
        "ERROR: Please assign the source programme with --prg. Possible "
        . "Options are 'exonerate', 'spaln', 'gth' or 'scipio'.\n";
    exit(1);
}

# check program source option
if (   $prgsrc ne "exonerate"
    && $prgsrc ne "spaln"
    && $prgsrc ne "gth"
    && $prgsrc ne "scipio"
    && $prgsrc ne "gemoma" )
{
    print STDERR
        "ERROR: Invalid value '$prgsrc' for option --prg. Possible Options "
        . "are 'exonerate', 'spaln', 'gth', 'scipio', or 'gemoma'.\n";
    exit(1);
}

# set source entry
if ( $prgsrc eq "exonerate" ) {
    $prgsrc = "xnt2h";
}

if ( $prgsrc eq "spaln" ) {
    $prgsrc           = "spn2h";
    $intron_threshold = $intron_threshold_spn;
}

if ( $prgsrc eq "gth" ) {
    $prgsrc           = "gth2h";
    $intron_threshold = $intron_threshold_gth;
}

if ( $prgsrc eq "scipio" ) {
    $prgsrc = "scipio2h";
}

if ( $prgsrc eq "gemoma" ) {
    $prgsrc = "gemoma2h";
}

if ( not( ( $prgsrc eq "xnt2h" ) || ( $prgsrc eq "gemoma2h" ) )
    && defined($genome_file) )
{
    print STDERR
        "ERROR: program name is $prgsrc and a genome file was specified. "
        . "Will ignore genome file.\n";
}
elsif ( $prgsrc eq "xnt2h" && defined($genome_file) ) {
    open( GENOME, "<", $genome_file )
        or die("Could not open genome fasta file $genome_file!\n");
    my $header;
    while (<GENOME>) {
        chomp;
        if (m/^>(.*)/) {
            $genome{$1} = "";
            $header = $1;
        }
        else {
            $genome{$header} .= $_;
        }
    }
    close(GENOME) or die("Could not close genome fasta file $genome_file!\n");
}

if ($CDS) {
    $CDSpartid = "CDS";
}

open( ALN,   "<$alignfile" )     or die("Cannot open file: $alignfile\n");
open( HINTS, ">$hintsfilename" ) or die("Cannot open file: $hintsfilename");

while (<ALN>) {
    chomp;
    my @f = split( /\t/, $_ );
    next unless ( scalar(@f) >= 8 );
    my $seqname = $f[0];
    my $type    = $f[2];
    my $start   = $f[3];
    my $end     = $f[4];
    my $score   = $f[5];
    my $strand  = $f[6];
    my $frame   = $f[7];
    if ( $end < $start ) {
        my $tmp = $start;
        $start = $end;
        $end   = $tmp;
    }

    # get target protein for exonerate
    if ( $type eq "gene" && $prgsrc eq "xnt2h" ) {
        /sequence (\S+) ; /;
        $parent = $1;
    }

    # get target protein for gth
    if ( $type eq "mRNA" && $prgsrc eq "gth2h" ) {
        my @info   = split( /\=/, $f[8] );
        my @rnaid  = split( /;/,  $info[1] );    # $rnaid[0]
        my @geneid = split( /;/,  $info[2] );    # $geneid[0]
        @info = split( /\s/, $info[-1] );
        $parent
            = $info[0] . "_" . $seqname . "_" . $rnaid[0] . "_" . $geneid[0];
    }

    # get target protein for spaln
    if ( $prgsrc eq "spn2h" && ( $type eq "CDS" || $type eq "cds" ) ) {
        my @info = split( /\=/, $f[8] );
        @info = split( /\s/, $info[-1] );
        $parent = $info[0];
    }

    # get target protein for scipio
    if ( $type eq "protein_match" && $prgsrc eq "scipio2h" ) {
        $f[8] =~ m/Query=(\S+) (\d+) (\d+)/;
        $parent = $1;
        $qstart = $2;    # start of alignment in query protein
        $qend   = $3;    # start of alignment in query protein
    }

    # get target protein for gemoma
    if ( $type eq "prediction" && $prgsrc eq "gemoma2h" ) {
        my @info = split( /;/, $f[8] );
        @info = split( /=/, $info[0] );
        $parent = $info[1];
    }

    # create start and stop hints from exonerate
    if ( ( $type eq "gene" && $prgsrc eq "xnt2h" && defined($genome_file) ) )
    { # || ($type eq "gene" && $prgsrc eq "gemoma2h" && defined($genome_file))){
        my $pot_start;
        if ( $strand eq "+" ) {
            $pot_start = substr( $genome{$seqname}, $start - 1, 3 );
        }
        elsif ( $strand eq "-" ) {
            $pot_start = substr( $genome{$seqname}, $end - 3, 3 );
            $pot_start =~ tr/acgtACGT/tgcaTGCA/;
            $pot_start = reverse($pot_start);
        }
        if ( defined($pot_start) ) {
            if ( $pot_start =~ m/(ATG)|(TTG)|(GTG)|(CTG)/i ) {
                print_start( $seqname, $strand, $start, $end );
            }
        }
        my $pot_stop;
        if ( $strand eq "+" ) {
            $pot_stop = substr( $genome{$seqname}, $end - 3, 3 );
        }
        else {
            $pot_stop = substr( $genome{$seqname}, $start - 1, 3 );
            $pot_stop =~ tr/acgtACGT/tgcaTGCA/;
            $pot_stop = reverse($pot_stop);
        }
        if ( $pot_stop =~ m/(TAA)|(TGA)|(TAG)/i ) {
            print_stop( $seqname, $strand, $start, $end );
        }
    }
    elsif ( $prgsrc eq "spn2h" || $prgsrc eq "gth2h" ) {
        if ((      $type eq "cds"
                && $f[8] =~ m/Target=.* (\d+) \d+ (\+|-)/
                && $prgsrc eq "spn2h"
            )
            || (   $type eq "mRNA"
                && $f[8] =~ m/Target=.* (\d+) \d+ (\+|-)/
                && $prgsrc eq "gth2h" )
            )
        {
            if ( $1 == 1 ) {
                print_start( $seqname, $strand, $start, $end );
                print_stop( $seqname, $strand, $start, $end );
            }
        }
    }
    elsif ( $prgsrc eq "gemoma2h" && $f[8] =~ m/start=(\w);stop=(.);/ ) {
        if ( $1 eq "M" ) {
            print_start( $seqname, $strand, $start, $end );
        }
        if ( $2 eq "*" ) {
            print_stop( $seqname, $strand, $start, $end );
        }
    }

    if ( $type eq "CDS" || $type eq "cds" || $type eq "protein_match" ) {
        if ( !$CDS ) {

            # CDSpart hint
            $start += $CDSpart_cutoff;
            $end -= $CDSpart_cutoff;
        }
        if ( $start > $end ) {
            $start = $end = int( ( $start + $end ) / 2 );
        }
        print HINTS
            "$seqname\t$prgsrc\t$CDSpartid\t$start\t$end\t$score\t$strand\t"
            . "$frame\tsrc=$source;grp=$parent;pri=$priority\n";
        if ( $prgsrc eq "spn2h" || $prgsrc eq "scipio2h" ) {
            get_intron( \@f );
        }
        elsif ( $prgsrc eq "gemoma2h" ) {
            get_gemoma_intron( \@f );
        }
    }

    if ( $type eq "exon" && $prgsrc eq "gth2h" ) {
        get_intron( \@f );
    }

    # intron hints for exonerate
    if ( $type eq "intron" ) {
        if (   $end - $start + 1 >= $minintronlen
            && $end - $start + 1 <= $maxintronlen )
        {
            print HINTS
                "$seqname\t$prgsrc\tintron\t$start\t$end\t$score\t$strand\t"
                . ".\tsrc=$source;grp=$parent;pri=$priority\n";
        }
    }
}
close(ALN)   or die("Could not close file $alignfile!\n");
close(HINTS) or die("Could not close file $hintsfilename!\n");

# intron hints for spaln and gth
sub get_intron {
    my $line = shift;
    $intron_score = $prevScore + @{$line}[5] / 2;
    if ( $prevParent ne $parent ) {
        if ( @{$line}[6] eq "-"
            && ( $prgsrc eq "spn2h" || $prgsrc eq "scipio2h" ) )
        {
            $intron_end = @{$line}[3] - 1
                ;   # these spliced aligners output in reverse order in genome
        }
        else {
            $intron_start = @{$line}[4] + 1;
        }
    }
    else {
        if ( @{$line}[6] eq "-"
            && ( $prgsrc eq "spn2h" || $prgsrc eq "scipio2h" ) )
        {
            $intron_start = @{$line}[4] + 1;
        }
        else {
            $intron_end = @{$line}[3] - 1;
        }
        if ( $intron_end < $intron_start ) {
            my $tmp = $intron_start;
            $intron_start = $intron_end;
            $intron_end   = $tmp;
        }

# check conditions: length of intron is at least $minintronlen and maximal
# $maxintronlen and its score is greater than $intron_threshold
        if (   $intron_end - $intron_start + 1 >= $minintronlen
            && $intron_end - $intron_start + 1 <= $maxintronlen
            && ( !defined($intron_threshold)
                || $intron_score > $intron_threshold )
            )
        {
            if ( $prgsrc ne "scipio2h"
                || ( defined($prevQend) && $prevQend + 1 == $qstart ) )
            {
                print HINTS
                    "@{$line}[0]\t$prgsrc\tintron\t$intron_start\t$intron_end"
                    . "\t$intron_score\t@{$line}[6]\t.\tsrc=$source;"
                    . "grp=$parent;pri=$priority\n";
            }
        }
        if ( @{$line}[6] eq "-"
            && ( $prgsrc eq "spn2h" || $prgsrc eq "scipio2h" ) )
        {
            $intron_end = @{$line}[3] - 1;
        }
        else {
            $intron_start = @{$line}[4] + 1;
        }
    }
    $prevScore  = @{$line}[5] / 2;
    $prevParent = $parent;
    $prevQend   = $qend if ( defined($qend) );
}

sub get_gemoma_intron {
    my $line = shift;
    if ( @{$line}[6] eq "+" ) {
        if ( $prevParent ne $parent ) {
            $intron_start = @{$line}[4] + 1;
        }
        else {
            $intron_end = @{$line}[3] - 1;
            if ( $intron_end < $intron_start ) {
                my $tmp = $intron_start;
                $intron_start = $intron_end;
                $intron_end   = $tmp;
            }
            if (   $intron_end - $intron_start + 1 >= $minintronlen
                && $intron_end - $intron_start + 1 <= $maxintronlen )
            {
                print HINTS
                    "@{$line}[0]\t$prgsrc\tintron\t$intron_start\t$intron_end"
                    . "\t.\t@{$line}[6]\t.\tsrc=$source;"
                    . "grp=$parent;pri=$priority\n";
            }
            $intron_start = @{$line}[4] + 1;
        }
    }
    else {
        if ( $prevParent ne $parent || not( defined($prevParent) ) ) {
            $intron_end = @{$line}[3] - 1;
        }
        else {
            $intron_start = @{$line}[4] + 1;
            if ( $intron_start > $intron_end ) {
                my $tmp = $intron_end;
                $intron_end   = $intron_start;
                $intron_start = $tmp;
            }
            if (   $intron_end - $intron_start + 1 >= $minintronlen
                && $intron_end - $intron_start + 1 <= $maxintronlen )
            {
                print HINTS
                    "@{$line}[0]\t$prgsrc\tintron\t$intron_start\t$intron_end"
                    . "\t.\t@{$line}[6]\t.\tsrc=$source;grp=$parent;"
                    . "pri=$priority\n";
            }
            $intron_end = @{$line}[3] - 1;
        }
    }
    $prevParent = $parent;
}

sub print_start {
    my $seqname = shift;
    my $strand  = shift;
    my $start   = shift;
    my $end     = shift;
    print HINTS "$seqname\t$prgsrc\tstart\t";
    if ( $strand eq "+" ) {
        print HINTS "$start\t" . ( $start + 2 );
    }
    else {
        print HINTS ( $end - 2 ) . "\t$end";
    }
    print HINTS "\t.\t$strand\t0\tsrc=$source;grp=$parent;pri=$priority\n";
}

sub print_stop {
    my $seqname = shift;
    my $strand  = shift;
    my $start   = shift;
    my $end     = shift;
    print HINTS"$seqname\t$prgsrc\tstop\t";
    if ( $strand eq "+" ) {
        print HINTS ( $end - 2 ) . "\t$end";
    }
    else {
        print HINTS "$start\t" . ( $start + 2 );
    }
    print HINTS "\t.\t$strand\t0\tsrc=$source;grp=$parent;pri=$priority\n";
}
