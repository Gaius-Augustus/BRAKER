#!/usr/bin/perl

####################################################################################################
#                                                                                                  #
# braker.pl                                                                                        #
# Pipeline for predicting genes with GeneMark-ET and AUGUSTUS with RNA-Seq                         #
#                                                                                                  #
# Authors: Simone Lange, Katharina Hoff, Mario Stanke                                              #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# Release date: January 7th 2015                                                                   #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
#                                                                                                  #
# Usage:                                                                                           #
# braker.pl [OPTIONS] --genome=genome.fa --bam=rnaseq.bam                                          #
#                                                                                                  #
####################################################################################################

# ----------------------------------------------------------------------
# | sub check_upfront              | autoAug.pl           | 07.01.2015 |
# | sub check_fasta_headers        | autoAug.pl           | 07.02.2015 |
# | helpMod qw(find chec...)       | helpMod.pm           | ??.??.???? |
# | first outline for braker       | Simone Lange         | 05.09.2014 |
# | uptodate integrated            |                      | 10.09.2014 |
# | print stdout,LOG output        |                      | 10.09.2014 |
# | .pm check for GeneMark         |                      | 10.09.2014 |
# | add parts from                 |                      | 30.09.2014 |
# | simplifyFastaHeaders.pl        | Katharina Hoff       | 03.12.2012 |
# | alteration on parts from       | Simone Lange         | 30.09.2014 |
# | simplifyFastaHeaders.pl &      |                      |            |
# | sub check_fasta_headers        |                      |            |
# | add filterIntronsFindStrand.pl |                      | 07.10.2014 |
# | add check whether augustus and | optimize_augustus.pl | 05.11.2014 |
# | etraining are executable       | (Mario Stanke)       | 23.04.2007 | 
# | add --optCfgfile, --fungus     | Simone Lange         | 10.11.2014 |
# | option, PATH also as variable  |                      |            |
# | fork AUGUSTUS prediction       |                      |            |
# | add parts from autoAugTrain.pl | Simone Lange         | 12.12.2014 |
# | parts from sub train           | autoAugTrain.pl      | 07.02.2015 |
# | TODO: add more options         |                      |            |
# ----------------------------------------------------------------------

use Getopt::Long;
use File::Compare;
use File::Path qw(make_path rmtree);
use Module::Load::Conditional qw(can_load check_install requires);
use Scalar::Util::Numeric qw(isint);
use POSIX qw(floor);
use List::Util qw(min);
use Parallel::ForkManager;

use Cwd;

use File::Spec::Functions qw(rel2abs); # $abs_path = File::Spec->rel2abs($path)
use File::Basename qw(dirname basename); # ?

BEGIN{
  $0 = rel2abs($0);
  our $directory = dirname($0);
} 
use lib $directory;
use helpMod qw(find checkFile formatDetector relToAbs setParInConfig uptodate);
use Term::ANSIColor qw(:constants); #?

use strict;
use warnings;


my $usage = <<'ENDUSAGE';

braker.pl     annotate genomes based on RNAseq using GeneMark-ET and AUGUSTUS

SYNOPSIS

braker.pl [OPTIONS] --genome=genome.fa --bam=rnaseq.bam


  --genome=genome.fa          fasta file with DNA sequences
  --bam=rnaseq.bam            bam file with spliced alignments from RNA-Seq


    
    
OPTIONS

    --help                               Print this help message
    --alternatives-from-evidence=true    Output alternative transcripts based on explicit evidence from 
                                         hints (default is true).
    --AUGUSTUS_CONFIG_PATH=/path/        Set path to AUGUSTUS (if not specified as environment variable).
      to/augustus                        Has higher priority than environment variable.
    --CPU                                Specifies the maximum number of CPUs that can be used during 
                                         computation
    --fungus                             GeneMark-ET option: run algorithm with branch point model (most 
                                         useful for fungal genomes)
    --GENEMARK_PATH=/path/to/            Set path to GeneMark-ET (if not specified as environment 
      gmes_petap.pl                      variable).
                                         Has higher priority than environment variable.
    --hints=hints.gff                    Alternatively to calling braker.pl with a bam file, it is 
                                         possible to call it with a file that contains introns extracted 
                                         from RNA-Seq data in gff format. This flag also allows the usage
                                         of hints from additional extrinsic sources for gene prediction 
                                         with AUGUSTUS. To consider such additional extrinsic information,
                                         you need to use the flag --optCfgFile to specify parameters for 
                                         all sources in the hints file
                                         (including the source "E" for intron hints from RNA-Seq).
    --optCfgFile=ppx.cfg                 Optional custom config file for AUGUSTUS (see --hints).
    --overwrite                          Overwrite existing files (except for species parameter files)
    --skipGeneMark-ET                    Skip GeneMark-ET and use provided GeneMark-ET output (e.g. from a
                                         different source) 
    --skipOptimize                       Skip optimize parameter step (not recommended).
    --species=sname                      Species name. Existing species will not be overwritten. 
                                         Uses Sp_1 etc., if no species is assigned
    --useexisting                        Use the present config and parameter files if they exist for 
                                         'species'
    --UTR                                Predict untranslated regions. Default is off.
    --workingdir=/path/to/wd/            Set path to working directory. In the working directory results
                                         and temporary files are stored
    --version                            print version number of braker.pl
                           

DESCRIPTION
      
  Example:

    braker.pl [OPTIONS] --genome=genome.fa  --species=speciesname --bam=accepted_hits.bam

ENDUSAGE

my $version = 1.2; # braker.pl version number
my $alternatives_from_evidence = "true"; # output alternative transcripts based on explicit evidence from hints
my $augpath;
my $augustus_cfg_path;                # augustus config path, higher priority than $AUGUSTUS_CONFIG_PATH on system
my $AUGUSTUS_CONFIG_PATH = $ENV{'AUGUSTUS_CONFIG_PATH'};
my @bam;                              # bam file names
my $bool_species = "true";            # false, if $species contains forbidden words (e.g. chmod)
my $cmdString;                        # to store shell commands
my $CPU = 1;                          # number of CPUs that can be used
my $currentDir = cwd();               # working superdirectory where program is called from
my $errorfile;                        # stores current error file name
my $errorfilesDir;                    # directory for error files
my @files;                            # contains all files in $rootDir
my $flanking_DNA;                     # length of flanking DNA, default value is min{ave. gene length/2, 10000}
my @forbidden_words;                  # words/commands that are not allowed in species name (e.g. unlink)
my $fungus = 0;                       # option for GeneMark-ET
my $gb_good_size;                     # number of LOCUS entries in 'genbank.good.gb'                         
my $genbank;                          # genbank file name
my $genemarkDir;                      # directory for GeneMark-ET output
my $GENEMARK_PATH = $ENV{'GENEMARK_PATH'}; # path to 'gmes_petap.pl' script on system
my $GMET_path;                        # GeneMark-ET path, higher priority than $GENEMARK_PATH
my $genome;                           # name of sequence file
my $genome_length = 0;                # length of genome file
my $help;                             # print usage
my @hints;                            # input hints file names
my $hintsfile;                        # hints file later used by AUGUSTUS
my $limit = 10000000;                 # maximum for generic species Sp_$limit
my $logfile;                          # contains used shell commands
my $optCfgFile;                       # optinonal extrinsic config file for AUGUSTUS
my $otherfilesDir;                    # directory for other files besides GeneMark-ET output and parameter files
my $overwrite = 0;                    # overwrite existing files (except for species parameter files)
my $parameterDir;                     # directory of parameter files for species
my $perlCmdString;                    # stores perl commands
my $printVersion = 0;
my $scriptPath=dirname($0);           # path of directory where this script is located
my $skipGeneMarkET = 0;               # skip GeneMark-ET and use provided GeneMark-ET output (e.g. from a different source) 
my $skipoptimize = 0;                 # skip optimize parameter step
my $species;                          # species name
my $stdoutfile;                       # stores current standard output
my $string;                           # string for storing script path
my $testsize;                         # AUGUSTUS training parameter: number of genes in a file that is
                                      # used to measure accuracy during parameter estimation with
                                      # optimize_augustus.pl. Default: 1000. If there are less than 1000
                                      # genes, all genes are used to measure accuracy. Decreasing this
                                      # parameter speeds up braker.pl but decreases prediction accuracy.
                                      # At least 300 genes are required for training AUGUSTUS.
my $useexisting = 0;                  # start with and change existing config, parameter and result files
my $UTR = "off";
my $workDir;                          # in the working directory results and temporary files are stored
my $bam;
my $hints;

# list of forbidden words for species name
@forbidden_words = ("system", "exec", "passthru", "run", "fork", "qx", "backticks", "chmod", "chown", "chroot", "unlink", "do", "eval", "kill", "rm", "mv", "grep", "cd", "top", "cp", "for", "done", "passwd", "while"); 


if(@ARGV==0){
  print "$usage\n"; 
  exit(0);
}

GetOptions( 'alternatives-from-evidence=s'  => \$alternatives_from_evidence,
            'AUGUSTUS_CONFIG_PATH=s'        => \$augustus_cfg_path,
            'bam=s'                         => \@bam,
            'CPU=i'                         => \$CPU,
            'fungus!'                       => \$fungus,
            'GENEMARK_PATH=s'               => \$GMET_path,
            'genome=s'                      => \$genome,
            'hints=s'                       => \@hints,
            'optCfgFile=s'                  => \$optCfgFile,
            'overwrite!'                    => \$overwrite,
            'skipGeneMark-ET!'              => \$skipGeneMarkET,
            'skipOptimize!'                 => \$skipoptimize,
            'species=s'                     => \$species,
            'testsize=i'                    => \$testsize,
            'useexisting!'                  => \$useexisting,
            'UTR=s'                         => \$UTR,
            'workingdir=s'                  => \$workDir,
            'help!'                         => \$help,
            'version!'                      => \$printVersion);

if($help){
  print $usage;
  exit(0);
}

if($printVersion){
    print "braker.pl version $version\n";
    exit(0);
}

             ############ make some regular checks ##############


# if no working directory is set, use current directory
if(!defined $workDir){
  $workDir = $currentDir;
}else{
  my $last_char = substr($workDir, -1);
  if($last_char eq "\/"){
    chop($workDir);
  }
}


# check the write permission of $workDir before building of the work directory
if(! -w $workDir){
  print STDERR "ERROR: Do not have write permission for $workDir.\nPlease use command 'chmod' to reset permission or specify another working directory\n"; # see autoAug.pl
  exit(1);
}


# set path to augustus config folder
if(defined($augustus_cfg_path)){
 $AUGUSTUS_CONFIG_PATH = $augustus_cfg_path;
}


# set path to GeneMark-ETs gmes_petap.pl script
if(defined($GMET_path)){
  $GENEMARK_PATH = $GMET_path;
}


# check upfront whether any common problems will occur later # see autoAug.pl
print STDOUT "NEXT STEP: check files and settings\n"; 
check_upfront();
check_options();


# check whether RNAseq files are specified
if(!@bam && !@hints){
  print STDERR "ERROR: No RNAseq or hints files specified. Please set at least one RNAseq file.\n$usage";
  exit(1);
}


# check whether bam files exists
if(@bam){
  @bam = split(/[\s,]/, join(',',@bam));
  for(my $i=0; $i<scalar(@bam); $i++){
    if(! -e $bam[$i]){
      print STDERR "ERROR: BAM file $bam[$i] does not exist. Please check.\n";
      exit(1);
    }
    $bam[$i] = rel2abs($bam[$i]);
  }
}


# check whether hints files exists
if(@hints){
  @hints = split(/[\s,]/, join(',',@hints));
  for(my $i=0; $i<scalar(@hints); $i++){
    if(! -e "$hints[$i]"){
      print STDERR "ERROR: Hints file $hints[$i] does not exist. Please check.\n";
      exit(1);
    }
    $hints[$i] = rel2abs($hints[$i]);
    check_gff($hints[$i]);
  }
}


# check whether species is specified
if(defined($species)){
  if($species =~ /[\s]/){
    print STDOUT "WARNING: Species name contains invalid white space characters. Will replace white spaces with underline character '_' \n";
    $species =~ s/\s/\_/g; print "species $species\n";
  }

  foreach my $word (@forbidden_words){
    if($species =~m/\A$word}\Z/){
      print STDOUT "WARNING: $species is not allowed as a species name. ";
      $bool_species = "false";
    }
  }
}
if(!defined($species) || $bool_species eq "false"){
  my $no = 1;
  $species = "Sp_$no";
  while($no <= $limit){
    $species = "Sp_$no";
    if((! -d "$AUGUSTUS_CONFIG_PATH/species/$species")){
      last;
    }else{
      $no++;
    }
  }
  if($bool_species eq "false"){
    print STDOUT "Program will use $species instead.\n";
  }else{
    print STDOUT "No species was set. Program will use $species.\n";
  }
}
    

# check species directory
if(-d "$AUGUSTUS_CONFIG_PATH/species/$species" && !$useexisting){
  print STDERR "ERROR: $AUGUSTUS_CONFIG_PATH/species/$species already exists. Choose another species name, delete this directory or use the existing species with the option --useexisting.\n";
  exit(1);
}

if(! -d "$AUGUSTUS_CONFIG_PATH/species/$species" && $useexisting){
  print STDOUT "WARNING: $AUGUSTUS_CONFIG_PATH/species/$species does not exist. Braker will create  the necessary files for species $species.\n";
}

# check whether $rootDir already exists
my $rootDir = "$workDir/braker";
if (-d "$rootDir/$species" && !$overwrite){
  print STDOUT "WARNING: $rootDir/$species already exists. Braker will use existing files, if they are newer than the input files. You can choose another working directory with --workingdir=dir or overwrite it with --overwrite\n";
}


# check whether genome file is set
if(!defined($genome)){
  print STDERR "ERROR: No genome file was specified. Please set a genome file.\n";
  exit(1);
}


# check whether genome file exist
if(! -f "$genome"){
  print STDERR "ERROR: Genome file $genome does not exist. Please check.\n";
  exit(1);
}else{
  # create $rootDir
  my $bool_rootDir = "false";
  if(! -d $rootDir){
    make_path($rootDir);
    $bool_rootDir = "true";
  }

  # set other directories
  $genemarkDir = "$rootDir/$species/GeneMark-ET";
  $parameterDir = "$rootDir/$species/species";
  $otherfilesDir = "$rootDir/$species";
  $errorfilesDir = "$rootDir/$species/errors";

  $logfile = "$otherfilesDir/braker.log";
  # create other directories if necessary
  my $bool_otherDir = "false";
  if(! -d $otherfilesDir){
    make_path($otherfilesDir);
    $bool_otherDir = "true";

  }
  open (LOG, ">>".$logfile) or die "Cannot open file: $logfile\n";
  if($bool_rootDir eq "true"){
    print LOG "\# ".localtime.": create working directory $rootDir\n";
    print LOG "mkdir $rootDir\n\n";
  }
  if($bool_otherDir eq "true"){
    print LOG "\# ".localtime.": create working directory $bool_otherDir\n";
    print LOG "mkdir $otherfilesDir\n\n";
  }
  if(! -d $genemarkDir){
    make_path($genemarkDir);
    print LOG "\# ".localtime.": create working directory $genemarkDir\n";
    print LOG "mkdir $genemarkDir\n\n";
  }
  # check whether genemark.gtf file exists, if skipGeneMark-ET option is used
  if($skipGeneMarkET){
    print STDOUT "REMARK: The GeneMark-ET step will be skipped.\n";
    if(! -f "$genemarkDir/genemark.gtf"){
      print STDERR "ERROR: The --skipGeneMark-ET option was used, but there is no genemark.gtf file under $genemarkDir. Please check or please do not skip the GeneMark-ET step.\n";
      exit(1);
    }
  }

  if(! -d $parameterDir){
    make_path($parameterDir);
    print LOG "\# ".localtime.": create working directory $parameterDir\n";
    print LOG "mkdir mkdir $parameterDir\n\n";
  }
  if(! -d $errorfilesDir){
    make_path($errorfilesDir);
    print LOG "\# ".localtime.": create working directory $errorfilesDir\n";
    print LOG "mkdir $errorfilesDir\n\n";
  }

  $genome = rel2abs($genome);
  $cmdString = "cd $rootDir;";
  print LOG "\# ".localtime.": change to working directory $rootDir\n";
  print LOG "$cmdString\n\n";
  chdir $rootDir or die ("Could not change to directory $rootDir.\n");

  check_fasta_headers($genome);
  make_hints();
  GeneMark_ET(); 
  new_species(); 
  training(); 
  augustus();
  clean_up();
  close(LOG) or die("Could not close log file $logfile!\n");
}



                           ############### sub functions ##############


         ####################### make hints #########################
# make hints from BAM files and optionally combine it with additional hints file
sub make_hints{
  my $bam_hints;
  my $hintsfile_temp = "$otherfilesDir/hintsfile.temp.gff";
  $hintsfile = "$otherfilesDir/hintsfile.gff";

  if(@bam){
    my $bam_temp = "$otherfilesDir/bam2hints.temp.gff";
    for(my $i=0; $i<scalar(@bam); $i++){
      $errorfile = "$errorfilesDir/bam2hints.$i.stderr";
      if(!uptodate([$bam[$i]],[$hintsfile])  || $overwrite){
        print STDOUT "NEXT STEP: make hints from BAM file $bam[$i]\n";
        $augpath = "$AUGUSTUS_CONFIG_PATH/../auxprogs/bam2hints/bam2hints";
        $cmdString = "$augpath --intronsonly --in=$bam[$i] --out=$bam_temp 2>$errorfile";
        print LOG "\# ".localtime.": make hints from BAM file $bam[$i]\n";
        print LOG "$cmdString\n\n";
        system("$cmdString")==0 or die("failed to execute: $!\n");
        $cmdString = "cat $bam_temp >>$hintsfile_temp";
        print LOG "\# ".localtime.": add hints from BAM file $bam[$i] to hints file\n";
        print LOG "$cmdString\n\n";
        system("$cmdString")==0 or die("failed to execute: $!\n");
        print STDOUT "hints from BAM file $bam[$i] added.\n";
      }
    }
    unlink($bam_temp);
  }

  if(@hints){
    for(my $i=0; $i<scalar(@hints); $i++){
      if(!uptodate([$hints[$i]],[$hintsfile]) || $overwrite){
        print STDOUT "NEXT STEP: add hints from file $hints[$i]\n";
        $cmdString = "cat $hints[$i] >> $hintsfile_temp";
        print LOG "\# ".localtime.": add hints from file $hints[$i]\n";
        print LOG "$cmdString\n\n";
        system("$cmdString")==0 or die("failed to execute: $!\n");
      }
    }
  }
  if(-f $hintsfile_temp || $overwrite){
    if(!uptodate([$hintsfile_temp],[$hintsfile]) || $overwrite){
      my $hintsfile_temp_sort = "$otherfilesDir/hints.temp.sort.gff";
      print STDOUT "NEXT STEP: sort hints\n";
      $cmdString = "cat $hintsfile_temp | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -n -k 3,3 | sort -s -k 1,1 >$hintsfile_temp_sort";
      print LOG "\# ".localtime.": sort hints\n";
      print LOG "$cmdString\n\n";
      system("$cmdString")==0 or die("failed to execute: $!\n");
      print STDOUT "hints sorted.\n";

      print STDOUT "NEXT STEP: summarize multiple identical hints to one\n";
      $string = find("join_mult_hints.pl");
      $errorfile = "$errorfilesDir/join_mult_hints.stderr";
      $perlCmdString = "perl $string <$hintsfile_temp_sort >$hintsfile_temp 2>$errorfile";
      print LOG "\# ".localtime.": join multiple hints\n";
      print LOG "$perlCmdString\n\n";
      system("$perlCmdString")==0 or die("failed to execute: $!\n");
      print STDOUT "hints joined.\n";
      unlink($hintsfile_temp_sort);
    }

    if(!uptodate([$hintsfile_temp],[$hintsfile]) || $overwrite){
      print STDOUT "NEXT STEP: filter introns, find strand and change score to \'mult\' entry\n";
      $string = find("filterIntronsFindStrand.pl");
      $errorfile = "$errorfilesDir/filterIntronsFindStrand.stderr";
      $perlCmdString = "perl $string $genome $hintsfile_temp --score 1>$hintsfile 2>$errorfile";
      print LOG "\# ".localtime.": filter introns, find strand and change score to \'mult\' entry\n";
      print LOG "$perlCmdString\n\n";
      system("$perlCmdString")==0 or die("failed to execute: $!\n");
      print STDOUT "strands found and score changed.\n";
      unlink($hintsfile_temp);
      print STDOUT "hints file complete.\n";
    }
    if(-z $hintsfile){
      print STDERR "ERROR: The hints file is empty. Maybe the genome and the RNA-seq file do not belong together.\n";
      exit(1);
    } 
  }
}



         ####################### GeneMark-ET #########################
# start GeneMark-ET and convert its output to real gtf format
sub GeneMark_ET{
  if(!$skipGeneMarkET){
    if(!uptodate([$genome,$hintsfile],["$genemarkDir/genemark.gtf"])  || $overwrite){
      $cmdString = "cd $genemarkDir";
      print LOG "\# ".localtime.": change to GeneMark-ET directory $genemarkDir\n";
      print LOG "$cmdString\n\n";
      chdir $genemarkDir or die ("Could not change to directory $genemarkDir.\n");

      print STDOUT "NEXT STEP: execute GeneMark-ET\n"; 
      $string = "$GENEMARK_PATH/gmes_petap.pl";
      $errorfile = "$errorfilesDir/GeneMark-ET.stderr";
      $stdoutfile = "$otherfilesDir/GeneMark-ET.stdout";
      $perlCmdString = "perl $string --sequence=$genome --ET=$hintsfile --cores=$CPU";
      if($fungus){
        $perlCmdString .= " --fungus";
      }
      $perlCmdString .= " 1>$stdoutfile 2>$errorfile";
      print LOG "\# ".localtime.": execute GeneMark-ET\n";
      print LOG "$perlCmdString\n\n";
      system("$perlCmdString")==0 or die("failed to execute: $perlCmdString\n");
      print STDOUT "GeneMark-ET finished.\n";
    }
  }

  # convert GeneMark-ET output to gtf format with doublequotes

  if(!uptodate(["$genemarkDir/genemark.gtf"],["$genemarkDir/genemark.c.gtf"])  || $overwrite){
    print STDOUT "NEXT STEP: convert GeneMark-ET to real gtf format\n"; 

    $string=find("filterGenemark.pl");
    $errorfile = "$errorfilesDir/filterGenemark.stderr";
    $stdoutfile = "$otherfilesDir/filterGenemark.stdout";
    $perlCmdString="perl $string --genemark=$genemarkDir/genemark.gtf --introns=$hintsfile --output=$genemarkDir/genemark.c.gtf 1>$stdoutfile 2>$errorfile";
    print LOG "\# ".localtime.": convert GeneMark-ET output to real gtf format\n";
    print LOG "$perlCmdString\n\n";
    system("$perlCmdString")==0 or die("failed to execute: $!\n");
    print STDOUT "GeneMark-ET conversion.\n";
    $cmdString = "cd $rootDir";
    print LOG "\# ".localtime.": change to working directory $rootDir\n";
    print LOG "$cmdString\n\n";
    chdir $rootDir or die ("Could not change to directory $rootDir.\n");
  }
}



         ####################### create a new species #########################
# create a new species $species and extrinsic file from generic one 
sub new_species{
  $augpath = "$AUGUSTUS_CONFIG_PATH/species/$species";
  if((!uptodate([$augpath."/$species\_metapars.cfg"],[$augpath."/$species\_parameters.cfg", $augpath."/$species\_exon_probs.pbl"]) && !$useexisting) || ! -d "$AUGUSTUS_CONFIG_PATH/species/$species"){
    print STDOUT "NEXT STEP: create new species $species\n"; 
    $string=find("new_species.pl");
    $errorfile = "$errorfilesDir/new_species.stderr";
    $perlCmdString="perl $string --species=$species 2>$errorfile";
    print LOG "\# ".localtime.": create new species $species\n";
    print LOG "$perlCmdString\n\n";
    system("$perlCmdString")==0 or die("failed to execute: $!\n");
  }

  # create extrinsic file
  my $extrinsic = "$AUGUSTUS_CONFIG_PATH/extrinsic/extrinsic.M.RM.E.W.cfg";
  my $extrinsic_cp = "$AUGUSTUS_CONFIG_PATH/species/$species/extrinsic.$species.cfg";
  if(! -e $extrinsic_cp){
    print STDOUT "NEXT STEP: create extrinsic file: $extrinsic_cp\n";
    print LOG "\# ".localtime.": create extrinsic file\n";
    print LOG "cp $AUGUSTUS_CONFIG_PATH/extrinsic/extrinsic.M.RM.E.W.cfg $AUGUSTUS_CONFIG_PATH/species/$species/extrinsic.$species.cfg\n\n"; 
    print STDOUT "species $species created.\n";
    open (EXTRINSIC, $extrinsic) or die "Cannot open file: $extrinsic\n";
    open (OUT, ">".$extrinsic_cp) or die "Cannot open file: $extrinsic_cp\n";
    my $GENERAL = "false";
    while(<EXTRINSIC>){
      chomp;
      next if($GENERAL eq "true" && $_ !~ m /^\#/);
      if($GENERAL eq "true" && $_ =~ m /^\#/){
        $GENERAL = "false";
      }
      print OUT "$_\n";
      
      if($_ =~ m/^\[GENERAL\]/){
        $GENERAL = "true";
        print OUT "\n";
        print OUT "      start        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print OUT "       stop        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print OUT "        tss        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print OUT "        tts        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print OUT "        ass        1   1    1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print OUT "        dss        1   1    1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print OUT "   exonpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1.003\n";
        print OUT "       exon        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print OUT " intronpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print OUT "     intron      .25       .2  M    1  1e+100  RM  1     1    E 1   50    W 1    1\n";
        print OUT "    CDSpart        1   1    1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print OUT "        CDS        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print OUT "    UTRpart        1   1    1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print OUT "        UTR        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print OUT "     irpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print OUT "nonexonpart        1        1  M    1  1e+100  RM  1     1.15 E 1    1    W 1    1\n";
        print OUT "  genicpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";

        print LOG "\# ".localtime.": edit extrinsic file and add\n$_\n";
        print LOG "\n";
        print LOG "      start        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print LOG "       stop        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print LOG "        tss        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print LOG "        tts        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print LOG "        ass        1   1    1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print LOG "        dss        1   1    1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print LOG "   exonpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1.003\n";
        print LOG "       exon        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print LOG " intronpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print LOG "     intron      .25       .2  M    1  1e+100  RM  1     1    E 1   50    W 1    1\n";
        print LOG "    CDSpart        1   1    1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print LOG "        CDS        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print LOG "    UTRpart        1   1    1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print LOG "        UTR        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print LOG "     irpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        print LOG "nonexonpart        1        1  M    1  1e+100  RM  1     1.15 E 1    1    W 1    1\n";
        print LOG "  genicpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
      }
    }
    close(OUT) or die("Could not close extrinsic file $extrinsic_cp!\n");
    close(EXTRINSIC) or die("Could not close extrinsic file $extrinsic!\n");
    print STDOUT "extrinsic file created.\n";
  }
}



         ####################### train AUGUSTUS #########################
# create genbank file and train AUGUSTUS parameters for given species
sub training{
  # get average gene length for flanking region
  print STDOUT "NEXT STEP: get flanking DNA size\n"; 
  open (GENELENGTH, "<$genemarkDir/genemark.average_gene_length.out") or die "Cannot open file: $genemarkDir/genemark.average_gene_length.out\n";
  my $average_length = <GENELENGTH>;
  close(GENELENGTH) or die("Could not close file $genemarkDir/genemark.average_gene_length.out!\n");
  $flanking_DNA = min((floor($average_length/2), 10000));
  # create genbank file from fasta input and GeneMark-ET output
  $string=find("gff2gbSmallDNA.pl");
  $genbank = "$otherfilesDir/genbank.gb";
  if(!uptodate([$genome,"$genemarkDir/genemark.c.gtf"],[$genbank])  || $overwrite){
    print STDOUT "NEXT STEP: create genbank file\n";
    $errorfile = "$errorfilesDir/gff2gbSmallDNA.stderr";
    if(-z "$genemarkDir/genemark.c.gtf"){
      print STDERR "ERROR: The GeneMark-ET output file is empty. Please check.\n";
      exit(1);
    }
    $perlCmdString = "perl $string $genemarkDir/genemark.c.gtf $genome $flanking_DNA $genbank 2>$errorfile";
    print LOG "\# ".localtime.": create genbank file\n";
    print LOG "$perlCmdString\n\n";
    system("$perlCmdString")==0 or die("failed to execute: $!\n");
    print STDOUT "genbank file created.\n";  
  }

  # filter genbank file and only use the genes considered "good" (i.e. genes whose introns are represented in hints file) 
  if(!uptodate([$genbank, "$genemarkDir/genemark.f.good.gtf"],["$otherfilesDir/genbank.good.gb"])  || $overwrite){
    print STDOUT "NEXT STEP: filter genbank file\n";
    $string=find("filterGenesIn_mRNAname.pl");
    $errorfile = "$errorfilesDir/filterGenesIn_mRNAname.stderr";
    if(-z "$genemarkDir/genemark.f.good.gtf"){
      print STDERR "ERROR: The GeneMark-ET output file does not contain genes that are supported by the intron hints.\n";
      exit(1);
    }
    $perlCmdString = "perl $string $genemarkDir/genemark.f.good.gtf $genbank 1>$otherfilesDir/genbank.good.gb 2>$errorfile";
    print LOG "\# ".localtime.": filter genbank file\n";
    print LOG "$perlCmdString\n\n";
    system("$perlCmdString")==0 or die("failed to execute: $!\n");
    print STDOUT "genbank file filtered.\n";
  }

  # split into training und test set
  if(!uptodate(["$otherfilesDir/genbank.good.gb"],["$otherfilesDir/genbank.good.gb.test", "$otherfilesDir/genbank.good.gb.train"])  || $overwrite){
    print STDOUT "NEXT STEP: split genbank file into train and test file\n";
    $string = find("randomSplit.pl");
    $errorfile = "$errorfilesDir/randomSplit.stderr";
    $gb_good_size = `grep -c LOCUS $otherfilesDir/genbank.good.gb`;
    if($gb_good_size < 300){
      print STDOUT "WARNING: Number of good genes is low ($gb_good_size). Recomended are at least 300 genes\n";
    }
    if($gb_good_size > 1000){
      $testsize = 1000;
      $perlCmdString = "perl $string $otherfilesDir/genbank.good.gb $testsize 2>$errorfile";
      print LOG "\# ".localtime.": split genbank file into train and test file\n";
      print LOG "$perlCmdString\n\n";
      system("$perlCmdString")==0 or die("failed to execute: $!\n");
      print STDOUT "genbank file splitted.\n";
    }
  }

  if(!$useexisting){
    # train AUGUSTUS for the first time
    if(!uptodate(["$otherfilesDir/genbank.good.gb.train","$otherfilesDir/genbank.good.gb"],["$otherfilesDir/firstetraining.stdout"])){
      # set "stopCodonExcludedFromCDS" to true
      print STDOUT "NEXT STEP: Seting value of \"stopCodonExcludedFromCDS\" in $AUGUSTUS_CONFIG_PATH/species/$species/$species\_parameters.cfg to \"true\"\n"; # see autoAugTrain.pl
      print LOG "\# ".localtime.": Seting value of \"stopCodonExcludedFromCDS\" in $AUGUSTUS_CONFIG_PATH/species/$species/$species\_parameters.cfg to \"true\"\n";
      setParInConfig($AUGUSTUS_CONFIG_PATH."/species/$species/$species\_parameters.cfg", "stopCodonExcludedFromCDS", "true"); # see autoAugTrain.pl

      # first try with etraining
      print STDOUT "NEXT STEP: first etraining\n"; 
      $augpath = "$AUGUSTUS_CONFIG_PATH/../bin/etraining";
      $errorfile = "$errorfilesDir/firstetraining.stderr";
      $stdoutfile = "$otherfilesDir/firstetraining.stdout";
      $gb_good_size = `grep -c LOCUS $otherfilesDir/genbank.good.gb`;
      if($gb_good_size <= 1000){
        $cmdString = "$augpath --species=$species $otherfilesDir/genbank.good.gb 1>$stdoutfile 2>$errorfile";
      }else{
        $cmdString = "$augpath --species=$species $otherfilesDir/genbank.good.gb.train 1>$stdoutfile 2>$errorfile";
      }
      print LOG "\# ".localtime.": first etraining\n";
      print LOG "$cmdString\n\n";
      system("$cmdString")==0 or die("failed to execute: $!\n");
      print STDOUT "first training complete.\n";

      # set "stopCodonExcludedFromCDS" to false and run etraining again if necessary
      my $t_b_t = $gb_good_size - $testsize;      # see autoAugTrain.pl
      my $err_stopCodonExcludedFromCDS = `grep -c "exon doesn't end in stop codon" $errorfile`; # see autoAugTrain.pl
      my $err_rate =  $err_stopCodonExcludedFromCDS / $t_b_t;  # see autoAugTrain.pl
      print STDOUT "Error rate of missing stop codon is $err_rate\n";  # see autoAugTrain.pl
      print LOG "\# ".localtime."Error rate of missing stop codon is $err_rate\n"; # see autoAugTrain.pl
      if($err_rate >= 0.5){ # see autoAugTrain.pl
        print STDOUT "The appropriate value for \"stopCodonExcludedFromCDS\" seems to be \"false\".\n"; # see autoAugTrain.pl
        print STDOUT "next step: Seting value of \"stopCodonExcludedFromCDS\" in $AUGUSTUS_CONFIG_PATH/species/$species/$species\_parameters.cfg to \"false\"\n"; # see autoAugTrain.pl
        setParInConfig($AUGUSTUS_CONFIG_PATH."/species/$species/$species\_parameters.cfg", "stopCodonExcludedFromCDS", "false");  # see autoAugTrain.pl
        print STDOUT "NEXT STEP: Trying etraining again\n";
        print LOG "\# ".localtime.": Trying etraining again\n";
        print LOG "$cmdString\n\n";
        system("$cmdString")==0 or die("failed to execute: $!\n");
        print STDOUT "trying etraining again complete.\n";
      }


      # adjust the stop-codon frequency in species_parameters.cfg according to train.out
      print LOG "\# ".localtime.": adjust the stop-codon frequency in species_parameters.cfg according to $stdoutfile\n";
      my $freqOfTag;            # see autoAugTrain.pl
      my $freqOfTaa;            # see autoAugTrain.pl
      my $freqOfTga;            # see autoAugTrain.pl
      open(TRAIN, "$stdoutfile") or die ("Can not open file $stdoutfile\n"); # see autoAugTrain.pl
      while(<TRAIN>){                 # see autoAugTrain.pl
        if(/tag:\s*.*\((.*)\)/){      # see autoAugTrain.pl
          $freqOfTag = $1;            # see autoAugTrain.pl
        }elsif(/taa:\s*.*\((.*)\)/){  # see autoAugTrain.pl
          $freqOfTaa = $1;            # see autoAugTrain.pl
        }elsif(/tga:\s*.*\((.*)\)/){  # see autoAugTrain.pl    
          $freqOfTga = $1;            # see autoAugTrain.pl
        }
      }
      close(TRAIN) or die("Could not close gff file $stdoutfile!\n");
      print LOG "\# ".localtime."Setting frequency of stop codons to tag=$freqOfTag, taa=$freqOfTaa, tga=$freqOfTga.\n";
      print STDOUT "NEXT STEP: Setting frequency of stop codons to tag=$freqOfTag, taa=$freqOfTaa, tga=$freqOfTga.\n";
      setParInConfig($AUGUSTUS_CONFIG_PATH."/species/$species/$species\_parameters.cfg", "/Constant/amberprob", $freqOfTag);  # see autoAugTrain.pl
      setParInConfig($AUGUSTUS_CONFIG_PATH."/species/$species/$species\_parameters.cfg", "/Constant/ochreprob", $freqOfTaa);  # see autoAugTrain.pl
      setParInConfig($AUGUSTUS_CONFIG_PATH."/species/$species/$species\_parameters.cfg", "/Constant/opalprob", $freqOfTga);   # see autoAugTrain.pl
      print STDOUT "frequency adjusted\n";
    }

    # first test
    if(!uptodate(["$otherfilesDir/genbank.good.gb.test", "$otherfilesDir/genbank.good.gb"],["$otherfilesDir/firsttest.stdout"])  || $overwrite){
      print STDOUT "NEXT STEP: first test\n";
      $augpath = "$AUGUSTUS_CONFIG_PATH/../bin/augustus";
      $errorfile = "$errorfilesDir/firsttest.stderr";
      $stdoutfile = "$otherfilesDir/firsttest.stdout";
      $gb_good_size = `grep -c LOCUS $otherfilesDir/genbank.good.gb`;
      if($gb_good_size <= 1000){
        $cmdString = "$augpath --species=$species $otherfilesDir/genbank.good.gb 1>$stdoutfile 2>$errorfile";
      }else{
        $cmdString = "$augpath --species=$species $otherfilesDir/genbank.good.gb.test 1>$stdoutfile 2>$errorfile";
      }
      print LOG "\# ".localtime.": first test\n";
      print LOG "$cmdString\n\n";
      system("$cmdString")==0 or die("failed to execute: $!\n");
      print STDOUT "first test finished.\n";
    }

    # optimize parameters
    if(!$skipoptimize){
      if(!uptodate(["$otherfilesDir/genbank.good.gb.train","$otherfilesDir/genbank.good.gb.test", "$otherfilesDir/genbank.good.gb"],[$AUGUSTUS_CONFIG_PATH."/species/$species/$species\_exon_probs.pbl", $AUGUSTUS_CONFIG_PATH."/species/$species/$species\_parameters.cfg", $AUGUSTUS_CONFIG_PATH."/species/$species/$species\_weightmatrix.txt"])){
        print STDOUT "NEXT STEP: optimize AUGUSTUS parameter\n";
        $string=find("optimize_augustus.pl");
        $errorfile = "$errorfilesDir/optimize_augustus.stderr";
        $stdoutfile = "$otherfilesDir/optimize_augustus.stdout";
        $gb_good_size = `grep -c LOCUS $otherfilesDir/genbank.good.gb`;
        if($gb_good_size <= 1000){
          $perlCmdString = "perl $string --species=$species --cpus=$CPU $otherfilesDir/genbank.good.gb 1>$stdoutfile 2>$errorfile";
        }else{
          $perlCmdString = "perl $string --species=$species --onlytrain=$otherfilesDir/genbank.good.gb.train --cpus=$CPU $otherfilesDir/genbank.good.gb.test 1>$stdoutfile 2>$errorfile";
        }
        print LOG "\# ".localtime.": optimize AUGUSTUS parameter\n";
        print LOG "$perlCmdString\n\n";
        system("$perlCmdString")==0 or die("failed to execute: $!\n");
        print STDOUT "parameter optimized.\n";
      }
    }

    # train AUGUSTUS for the second time
    if(!uptodate(["$otherfilesDir/genbank.good.gb.train","$otherfilesDir/genbank.good.gb"],["$otherfilesDir/secondetraining.stdout"])){
      print STDOUT "NEXT STEP: second etraining\n";
      $augpath = "$AUGUSTUS_CONFIG_PATH/../bin/etraining";
      $errorfile = "$errorfilesDir/secondetraining.stderr";
      $stdoutfile = "$otherfilesDir/secondetraining.stdout";
      $gb_good_size = `grep -c LOCUS $otherfilesDir/genbank.good.gb`;
      if($gb_good_size <= 1000){
        $cmdString = "$augpath --species=$species $otherfilesDir/genbank.good.gb 1>$stdoutfile 2>$errorfile";
      }else{
        $cmdString = "$augpath --species=$species $otherfilesDir/genbank.good.gb.train 1>$stdoutfile 2>$errorfile";
      }
      print LOG "\# ".localtime.": second etraining\n";
      print LOG "$cmdString\n\n";
      system("$cmdString")==0 or die("failed to execute: $!\n");
      print STDOUT "second etraining complete\n";
    }

    # second test
    if(!uptodate(["$otherfilesDir/genbank.good.gb.test","$otherfilesDir/genbank.good.gb"],["$otherfilesDir/secondtest.out"]) || $overwrite){
      print STDOUT "NEXT STEP: second test\n";
      $errorfile = "$errorfilesDir/secondtest.stderr";
      $stdoutfile = "$otherfilesDir/secondtest.stdout";
      $gb_good_size = `grep -c LOCUS $otherfilesDir/genbank.good.gb`;
      if($gb_good_size <= 1000){
        $cmdString = "$augpath --species=$species $otherfilesDir/genbank.good.gb >$stdoutfile 2>$errorfile";
      }else{
        $cmdString = "$augpath --species=$species $otherfilesDir/genbank.good.gb.test >$stdoutfile 2>$errorfile";
      }
      print LOG "\# ".localtime.": second test\n";
      print LOG "$cmdString\n\n";
      system("$cmdString")==0 or die("failed to execute: $!\n");
      print STDOUT "second test finished.\n";  
    }
  }
  
  # copy species files to working directory
  if(! -d "$parameterDir/$species"){
    print STDOUT "NEXT STEP: copy optimized parameters to working directory\n";
    $cmdString = "cp -r $AUGUSTUS_CONFIG_PATH/species/$species $parameterDir";
    print LOG "\# ".localtime.": copy optimized parameters to working directory\n";
    print LOG "$cmdString\n\n";
    system("$cmdString")==0 or die("failed to execute: $!\n");
    print STDOUT "parameter files copied.\n";
  }
}



         ####################### run AUGUSTUS and convert its output #########################
# run AUGUSTUS for given species with given options
sub augustus{
  $augpath = "$AUGUSTUS_CONFIG_PATH/../bin/augustus";
  my $scriptpath = "$AUGUSTUS_CONFIG_PATH/../scripts";
  my $extrinsic;
  my @genome_files;
  my $pm;
  if(!$useexisting){
    $extrinsic = "$AUGUSTUS_CONFIG_PATH/species/$species/extrinsic.$species.cfg";
  }else{
    $extrinsic = "$AUGUSTUS_CONFIG_PATH/extrinsic/extrinsic.M.RM.E.W.cfg";
  }
  if(!uptodate([$extrinsic,$hintsfile, $genome],["$otherfilesDir/augustus.gff"])  || $overwrite){
    if($CPU > 1){
      print STDOUT "NEXT STEP: split genome file in smaller parts\n";
      $string = find("splitMfasta.pl");    
      $errorfile = "$errorfilesDir/splitMfasta.stderr";
      my $minsize = floor($genome_length / $CPU);
      $perlCmdString = "perl $string $genome --outputpath=$otherfilesDir --minsize=$minsize 2>$errorfile";
      system("$perlCmdString")==0 or die("failed to execute: $!\n");
      @genome_files = `find $otherfilesDir -name "genome.split.*"`;
      print LOG "\# ".localtime.": split genome file in ".scalar(@genome_files)." parts\n";
      print LOG "$perlCmdString\n\n";
      print STDOUT "split genome file in ".scalar(@genome_files)." parts complete.\n";
    }else{
      push(@genome_files, $genome);
    }
    $pm = new Parallel::ForkManager($CPU);
    print STDOUT "NEXT STEP: run AUGUSTUS\n";
    for(my $i = 0; $i < scalar(@genome_files); $i++){
      chomp($genome_files[$i]);
      my $pid = $pm->start and next;
      my $idx = $i + 1;
      $errorfile = "$errorfilesDir/augustus.$idx.stderr";
      $stdoutfile = "$otherfilesDir/augustus.$idx.gff";
      $cmdString = "$augpath --species=$species --extrinsicCfgFile=$extrinsic --hintsfile=$hintsfile --UTR=$UTR";
      if(defined($optCfgFile)){
        $cmdString .= " --optCfgFile=$optCfgFile"; 
      }
      $cmdString .= " $genome_files[$i] 1>$stdoutfile 2>$errorfile";
      print LOG "\# ".localtime.": run AUGUSTUS for file $genome_files[$i]\n";
      print LOG "$cmdString\n\n";
      system("$cmdString")==0 or die("failed to execute: $!\n");
      $pm->finish;

    }
    $pm->wait_all_children;
    print STDOUT "AUGUSTUS complete.\n";

    if($CPU > 1){
      print STDOUT "NEXT STEP: concatenate and join AUGUSTUS output files\n";
      $string = find("join_aug_pred.pl");
      $cmdString = "cat $otherfilesDir/augustus.*gff | $string >$otherfilesDir/augustus.gff";
      print LOG "\# ".localtime.": concatenate and join AUGUSTUS output files\n";
      print LOG "$cmdString\n\n";
      system("$cmdString")==0 or die("failed to execute: $!\n");
      print STDOUT "join AUGUSTUS predictions complete.\n";
      for(my $idx = 1; $idx <= scalar(@genome_files); $idx++){
        unlink("$otherfilesDir/augustus.$idx.gff");
      }
    }else{
      $cmdString = "mv $otherfilesDir/augustus.1.gff $otherfilesDir/augustus.gff";
      print LOG "\# ".localtime.": rename AUGUSTUS file\n";
      print LOG "$cmdString\n\n";
      system("$cmdString")==0 or die("failed to execute: $!\n");
    }
  }
}

   
  

         ####################### delete all zero sized files #########################
# delete empty files
sub clean_up{
  print STDOUT "NEXT STEP: delete empty files\n";
  @files = `find $otherfilesDir -empty`;
  print LOG "\# ".localtime.": delete empty files\n";
  for(my $i=0; $i <= $#files; $i++){
    chomp($files[$i]); # to prevent error: Unsuccessful stat on filename containing newline
    if(-f $files[$i]){
      print LOG "rm $files[$i]\n";
      unlink(rel2abs($files[$i]));
    }
  }
  print STDOUT "empty files deleted.\n";
}



         ########################### some checks beforehand ############################
# check upfront whether any common problems will occur later
# find out if some programs are not installed.
# TODO: put more checks in here: blat, samtools, gmap, tophat 
# checks for GeneMark-ET: perl modules: YAML, Hash::Merge, Logger::Simple, Parallel::ForkManager
# checks for braker: perl modules: Scalar::Util::Numeric
sub check_upfront{ # see autoAug.pl
  my $pmodule;
  if(!$ENV{'AUGUSTUS_CONFIG_PATH'} && !defined($augustus_cfg_path)){ # see autoAug.pl
    print STDERR "ERROR: The environment variable AUGUSTUS_CONFIG_PATH is not defined. Please export an environment variable for AUGUSTUS or use --augustus-cfg-path=path/to/augustus.\n"; # see autoAug.pl
    exit(1);
  }
  
  if(!$ENV{'GENEMARK_PATH'}){
    print STDERR "ERROR: The environment variable GENEMARK_PATH to the 'gmes_petap.pl' script is not defined. Please export an environment variable or use --genemark-ET-path=path/to/gemes_petap.pl.\n"; 
    exit(1);
  }

  $augpath = "$AUGUSTUS_CONFIG_PATH/../bin/augustus";
  if(system("$augpath > /dev/null 2> /dev/null") != 0){                   # see autoAug.pl
    if(! -f $augpath){                                                    # see autoAug.pl
      print STDERR "ERROR: augustus executable not found at $augpath.\n"; # see autoAug.pl
    }else{
      print STDERR "ERROR: $augpath not executable on this machine.\n";   # see autoAug.pl
    }
    exit(1);
  }

  my $etrainpath = "$AUGUSTUS_CONFIG_PATH/../bin/etraining";
  if(system("$etrainpath > /dev/null 2> /dev/null") != 0){                   
    if(! -f $etrainpath){                                                    
      print STDERR "ERROR: etraining executable not found at $etrainpath.\n";
    }else{
      print STDERR "ERROR: $etrainpath not executable on this machine.\n";
    }
    exit(1);
  }

  find("gff2gbSmallDNA.pl");
  find("filterGenemark.pl");
  find("filterIntronsFindStrand.pl");
  find("new_species.pl");
  find("filterGenesIn_mRNAname.pl");
  find("join_mult_hints.pl");
  find("randomSplit.pl");
  find("optimize_augustus.pl");
  find("splitMfasta.pl");
  find("join_aug_pred.pl");

  # check whether required perl modules are installed (perl modules: YAML, Hash::Merge, Logger::Simple, Parallel::ForkManager)
  my $module_list = {
    "YAML",
    "Hash::Merge",
    "Logger::Simple",
    "Parallel::ForkManager",
  };

  $pmodule = check_install(module => "YAML");
  if(!$pmodule){
    print STDOUT "WARNING: Perl module 'YAML' is required but not installed yet.\n";
  }

  $pmodule = check_install(module => "Hash::Merge");
  if(!$pmodule){
    print STDOUT "WARNING: Perl module 'Hash::Merge' is required but not installed yet.\n";
  }

  $pmodule = check_install(module => "Logger::Simple");
  if(!$pmodule){
    print STDOUT "WARNING: Perl module 'Logger::Simple' is required but not installed yet.\n";
  }

  $pmodule = check_install(module => "Parallel::ForkManager");
  if(!$pmodule){
    print STDOUT "WARNING: Perl module 'Parallel::ForkManager' is required but not installed yet.\n";
  }

  $pmodule = check_install(module => "Scalar::Util::Numeric");
  if(!$pmodule){
    print STDOUT "WARNING: Perl module 'Scalar::Util::Numeric' is required but not installed yet.\n";
  }
}



# check whether hints file is in gff format
sub check_gff{
  my $gfffile = shift;
  print STDOUT "NEXT STEP: check if input file $gfffile is in gff format\n";
  open (GFF, $gfffile) or die "Cannot open file: $gfffile\n";
  while(<GFF>){
    my @gff_line = split(/\t/, $_);
    if(scalar(@gff_line) != 9){
      print STDERR "ERROR: File $gfffile is not in gff format. Please check\n";
      close(GFF) or die("Could not close gff file $gfffile!\n");
      exit(1);
    }else{
      if(!isint($gff_line[3]) || !isint($gff_line[4]) || $gff_line[5] =~ m/[^\d\.]/g || $gff_line[6] !~ m/[\+\-\.]/ || length($gff_line[6]) != 1 || $gff_line[7] !~ m/[0-2\.]{1}/ || length($gff_line[7]) != 1){
        print STDERR "ERROR: File $gfffile is not in gff format. Please check\n";
        close(GFF) or die("Could not close gff file $gfffile!\n");
        exit(1);
      }
    }
  }
  close(GFF) or die("Could not close gff file $gfffile!\n");
  print STDOUT "gff format check complete.\n";
}

# check whether all options are set correctly
sub check_options{
  print STDOUT "NEXT STEP: check options\n";
  if($alternatives_from_evidence ne "true" && $alternatives_from_evidence ne "false"){
    print STDERR "ERROR: \"$alternatives_from_evidence\" is not a valid option for --alternatives-from-evidence. Please use either true or false.\n";
    exit(1);
  } 

  if($UTR ne "on" && $UTR ne "off"){
    print STDERR "ERROR: \"$UTR\" is not a valid option for --UTR. Please use either on or off.\n";
    exit(1);
  }

  if(!isint($CPU)){
    print STDERR "ERROR: \"$CPU\" is not a valid option for --CPU. Please use an integer.\n";
    exit(1);
  }else{
    my $cpus_available = `nproc`;
    if($cpus_available < $CPU){
      print STDOUT "WARNING: Your system does not have $CPU cores available, only $cpus_available. Braker will use the $cpus_available available instead of the chosen $CPU.\n";
    }
  }
  print STDOUT "... options check complete.\n";
}

# check fasta headers
sub check_fasta_headers{            # see autoAug.pl
  
  my $fastaFile = shift;            # see autoAug.pl
  my $someThingWrongWithHeader = 0; # see autoAug.pl
  my $spaces = 0;                   # see autoAug.pl
  my $orSign = 0;                   # see autoAug.pl
  my $emptyC = 0;                   # see simplifyFastaHeaders.pl
  my $wrongNL = 0;                  # see simplifyFastaHeaders.pl
  my $prot = 0;                     # see simplifyFastaHeaders.pl
  my $dna = 0;                      # see simplifyFastaHeaders.pl
  my $mapFile = "$otherfilesDir/header.map";
  my $stdStr = "This may later on cause problems! The pipeline will create a new file without spaces and a header.map file to look up the old and new headers. This message will be suppressed from now on!\n";
  if(!uptodate([$genome],["$otherfilesDir/genome.fa"]) || $overwrite){
    print STDOUT "NEXT STEP: check fasta headers\n";
    open(FASTA, "<", $fastaFile) or die("Could not open fasta file $fastaFile!\n");
    open(OUTPUT, ">", "$otherfilesDir/genome.fa") or die("Could not open fasta file $otherfilesDir/genome.fa!\n");
    open(MAP, ">", $mapFile) or die("Could not open map file $mapFile.\n");
    while(<FASTA>){
      # check newline character
      if(not($_ =~ m/\n$/)){  # see simplifyFastaHeaders.pl
        if($wrongNL < 1){     # see simplifyFastaHeaders.pl
          print STDOUT "WARNING: something seems to be wrong with the newline character! This is likely to cause problems with the braker.pl pipeline and the AUGUSTUS web service! Please adapt your file to UTF8! This warning will be supressed from now on!\n"; # see simplifyFastaHeaders.pl
          $wrongNL++;         # see simplifyFastaHeaders.pl
        }
      }
      chomp;
      # look for whitespaces in fasta file
      if($_ =~ m/\s/){      # see autoAug.pl
        if($spaces == 0){   # see autoAug.pl
          print STDOUT "WARNING: Detected whitespace in fasta header of file $fastaFile. ".$stdStr; # see autoAug.pl
          $spaces++;        # see autoAug.pl
        }
      }

      # look for | in fasta file
      if($_ =~ m/\|/){      # see autoAug.pl
        if($orSign == 0){   # see autoAug.pl
          print STDOUT "WARNING: Detected | in fasta header of file $fastaFile. ".$stdStr; # see autoAug.pl
          $orSign++;        # see autoAug.pl
        }
      }

      # look for special characters in headers
      if(($_ !~ m/[>a-zA-Z0-9]/) && ($_ =~ m/^>/) ){
        if($someThingWrongWithHeader==0){   # see autoAug.pl
          print STDOUT "WARNING: Fasta headers in file $fastaFile seem to contain non-letter and non-number characters. That means they may contain some kind of special character. ".$stdStr; # see autoAug.pl
          $someThingWrongWithHeader++;      # see autoAug.pl
        }
      }

      if($_ =~ m/^>/){
        my $line = substr($_,1);
        # remove whitespaces, if necessary
        my @fasta_line = split(/\s/, $line);
        print OUTPUT ">$fasta_line[0]\n";   # see simplifyFastaHeaders.pl
        print MAP ">$fasta_line[0]\t$_\n";  # see simplifyFastaHeaders.pl  
      }else{
        if(length($_) > 0){    # see simplifyFastaHeaders.pl
          $genome_length += length($_);
          print OUTPUT "$_\n"; # see simplifyFastaHeaders.pl
          if($_ !~ m/[ATGCNatgcn]/){ # see simplifyFastaHeaders.pl
            if($dna == 0){ # see simplifyFastaHeaders.pl
              print STDOUT "Assuming that this is not a DNA fasta file because other characters than A, T, G, C, N, a, t, g, c, n were contained. If this is supposed to be a DNA fasta file, check the content of your file! If this is supposed to be a protein fasta file, please ignore this message!\n"; # see simplifyFastaHeaders.pl
              $dna++;      # see simplifyFastaHeaders.pl
            }
          }
          if($_ !~ m/[AaRrNnDdCcEeQqGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzJjXx]/){ # see simplifyFastaHeaders.pl
            if($prot == 0){ # see simplifyFastaHeaders.pl
              print STDOUT "Assuming that this is not a protein fasta file because other characters than AaRrNnDdCcEeQqGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzJjXx were contained. If this is supposed to be DNA fasta file, please ignore this message.\n"; # see simplifyFastaHeaders.pl
              $prot++;      # see simplifyFastaHeaders.pl
            }
          }
        }else{
          if($emptyC < 1){  # see simplifyFastaHeaders.pl
            print STDOUT "WARNING: empty line was removed! This warning will be supressed from now on!\n"; # see simplifyFastaHeaders.pl
          } 
          $emptyC++;        # see simplifyFastaHeaders.pl
        }
      }
    }
    close(FASTA) or die("Could not close fasta file $fastaFile!\n");
    close(OUTPUT) or die("Could not close output fasta file $otherfilesDir/genome.fa!\n");
    close(MAP) or die("Could not close map file $mapFile!\n");
    print STDOUT "fasta headers check complete.\n";
  } 
  $genome = "$otherfilesDir/genome.fa";
}

