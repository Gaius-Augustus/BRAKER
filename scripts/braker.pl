#!/usr/bin/env perl

####################################################################################################
#                                                                                                  #
# braker.pl                                                                                        #
# Pipeline for predicting genes with GeneMark-ET and AUGUSTUS with RNA-Seq                         #
#                                                                                                  #
# Authors: Simone Lange, Katharina Hoff, Mario Stanke                                              #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# Release date: April 26th 2017                                                                    #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
#                                                                                                  #
# Usage:                                                                                           #
# braker.pl [OPTIONS] --genome=genome.fa --bam=rnaseq.bam                                          #
#                                                                                                  #
####################################################################################################

use Getopt::Long;
use File::Compare;
use File::Path qw(make_path rmtree);
use Module::Load::Conditional qw(can_load check_install requires);
use Scalar::Util::Numeric qw(isint);
use POSIX qw(floor);
use List::Util qw(min);
use Parallel::ForkManager;

use Cwd;
use Cwd 'abs_path';

use File::Spec::Functions qw(rel2abs);
use File::Basename qw(dirname basename);

BEGIN{
  $0 = rel2abs($0);
  our $directory = dirname($0);
} 
use lib $directory;
use helpMod qw(find checkFile formatDetector relToAbs setParInConfig uptodate);
use Term::ANSIColor qw(:constants);

use strict;
use warnings;


my $usage = <<'ENDUSAGE';

braker.pl     Pipeline for predicting genes with GeneMark-ET and AUGUSTUS with RNA-Seq

SYNOPSIS

braker.pl [OPTIONS] --genome=genome.fa --bam=rnaseq.bam


  --genome=genome.fa          fasta file with DNA sequences
  --bam=rnaseq.bam            bam file with spliced alignments from RNA-Seq


    
    
OPTIONS

    --help                               Print this help message
    --nice                               Execute all system calls within braker.pl and its submodules
                                         with bash "nice" (default nice value)
    --alternatives-from-evidence=true    Output alternative transcripts based on explicit evidence from 
                                         hints (default is true).
    --augustus-args="--some_arg=bla"     One or several command line arguments to be passed to AUGUSTUS,
                                         if several arguments are given, separate by whitespace, i.e.
                                         "--first_arg=sth --second_arg=sth". 
    --AUGUSTUS_CONFIG_PATH=/path/        Set path to config directory of AUGUSTUS (if not specified as 
                                         environment variable). BRAKER1 will assume that the directories
                                         ../bin and ../scripts of AUGUSTUS are located relative to
                                         the AUGUSTUS_CONFIG_PATH. If this is not the case, please
                                         specify AUGUSTUS_BIN_PATH (and AUGUSTUS_SCRIPTS_PATH if required).
                                         The Perl commandline argument --AUGUSTUS_CONFIG_PATH has higher
                                         priority than the environment variable with the same name.
    --AUGUSTUS_BIN_PATH=/path/           Set path to the AUGUSTUS directory that contains binaries, i.e.
    					 augustus and etraining. This variable must only be set if 
    					 AUGUSTUS_CONFIG_PATH does not have ../bin and ../scripts of
    					 AUGUSTUS relative to its location i.e. for global AUGUSTUS 
    					 installations. BRAKER1 will assume that the directory
    					 ../scripts of AUGUSTUS is located relative to the AUGUSTUS_BIN_PATH.
    					 If this is not the case, please specify AUGUSTUS_SCRIPTS_PATH.
    --AUGUSTUS_SCRIPTS_PATH=/path/       Set path to AUGUSTUS directory that contains scripts, i.e.
    				         splitMfasta.pl. This variable most only be set if
    				         AUGUSTUS_CONFIG_PATH or AUGUSTUS_BIN_PATH do not contains the
    				         ../scripts directory of AUGUSTUS relative to their location, i.e.
    				         for special cases of a global AUGUSTUS installation.
    --BAMTOOLS_PATH=/path/to/            Set path to bamtools (if not specified as environment 
      bamtools/                          variable). Has higher priority than the environment variable.
    --cores                              Specifies the maximum number of cores that can be used during 
                                         computation
    --extrinsicCfgFile                   Optional. This file contains the list of used sources for the 
                                         hints and their boni and mali. If not specified the file "extrinsic.cfg" 
                                         in the config directory $AUGUSTUS_CONFIG_PATH is copied and adjusted.
    --fungus                             GeneMark-ET option: run algorithm with branch point model (most 
                                         useful for fungal genomes)
    --GENEMARK_PATH=/path/to/            Set path to GeneMark-ET (if not specified as environment 
      gmes_petap.pl/                     variable). Has higher priority than environment variable.
    --gff3                               Output in GFF3 format.
    --rnaseq-hints=hints.gff             Alternatively to calling braker.pl with a bam file, it is 
                                         possible to call it with a file that contains introns extracted 
                                         from RNA-Seq data in gff format. This flag also allows the usage
                                         of hints from additional extrinsic sources for gene prediction 
                                         with AUGUSTUS. To consider such additional extrinsic information,
                                         you need to use the flag --optCfgFile to specify parameters for 
                                         all sources in the hints file
                                         (including the source "E" for intron hints from RNA-Seq).
    --optCfgFile=ppx.cfg                 Optional custom config file for AUGUSTUS (see --rnaseq-hints).
    --overwrite                          Overwrite existing files (except for species parameter files)
    --SAMTOOLS_PATH=/path/to/            Optionally set path to samtools (if not specified as environment 
      samtools/                          variable) to fix BAM files automatically, if necessary. Has higher     
                                         priority than environment variable.
    --skipGeneMark-ET                    Skip GeneMark-ET and use provided GeneMark-ET output (e.g. from a
                                         different source) 
    --skipOptimize                       Skip optimize parameter step (not recommended).
    --softmasking                        Softmasking option for soft masked genome files. Set to 'on' or '1'
    --species=sname                      Species name. Existing species will not be overwritten. 
                                         Uses Sp_1 etc., if no species is assigned                          
    --useexisting                        Use the present config and parameter files if they exist for 
                                         'species'
    --UTR                                Predict untranslated regions. Default is off. (Only works if UTR
                                         parameters were previously optimized outside of BRAKER!)
    --workingdir=/path/to/wd/            Set path to working directory. In the working directory results
                                         and temporary files are stored
    --filterOutShort                     It may happen that a "good" training gene, i.e. one that has intron
                                         support from RNA-Seq in all introns predicted by GeneMark, is in fact
                                         too short. This flag will discard such genes that have supported introns
                                         and a neighboring RNA-Seq supported intron upstream of the start codon 
                                         within the range of the maximum CDS size of that gene and with a 
                                         multiplicity that is at least as high as 20% of the average intron 
                                         multiplicity of that gene.
    --crf                                Execute CRF training for AUGUSTUS; resulting parameters are only kept for
                                         final predictions if they show higher accuracy than HMM parameters.
    --prot_seq=prot.fa                   A protein sequence file in multiple fasta format. This file will be used
                                         to generate protein hints for AUGUSTUS by running one of the three
                                         alignment tools Exonerate (--prg=exonerate), Spaln (--prg=spaln) or 
                                         GenomeThreader (--prg=gth). Default is GenomeThreader if the tool is not 
                                         specified.
                                         Currently, hints from proteins are only used in the prediction step with
                                         AUGUSTUS. If --protein-seq is specified, it is not allowed to 
                                         specify --protein-aln or --protein-hints.
    --prot_aln=prot.aln                  Alignment file generated from aligning protein sequences against the 
                                         genome with either Exonerate (--prg=exonerate), or Spaln (--prg=spaln), or
                                         GenomeThreader (--prg=gth).
                                         To prepare alignment file, run Spaln2 with the following command:
                                            spaln -O0 ... > spalnfile
                                         To prepare alignment file, run Exonerate with the following command:
                                            exonerate --model protein2genome --showtargetgff T ... > exfile
                                         To prepare alignment file, run GenomeThreader with the following command:
                                            gth -genomic genome.fa  -protein protein.fa -gff3out -skipalignmentout 
                                               ... -o gthfile
                                         A valid option prg=... must be specified in combination with 
                                         --prot_aln. Generating tool will not be guessed.
                                         Currently, hints from proteins are only used in the prediction step with
                                         AUGUSTUS.
    --prot_hints   =prothints.gff        File containing protein hints for gene prediction with AUGUSTUS.
                                         Allowed features: NEED TO ADD THE FEATURES, LATER!
    --prg=gth|exonerate|spaln            Alignment tool ght (GenomeThreader), exonerate (Exonerate) or Spaln2
                                         (spaln) that will be used to generate protein alignments that will be the 
                                         basis for hints generation for gene prediction with AUGUSTUS (if specified
                                         in combination with --prot_seq) or that was used to externally
                                         generate an alignment file with the commands listed in description of 
                                         --prot_aln (if used in combination with --prot_aln).
    --ALIGNMENT_TOOL_PATH=/path/to/tool  Set path to alignment tool (GenomeThreader, Spaln, or Exonerate) if not 
                                         specified as environment variable. Has higher priority than environment
                                         variable.
    --version                            print version number of braker.pl
                           

DESCRIPTION
      
  Example:

    braker.pl [OPTIONS] --genome=genome.fa  --species=speciesname --bam=accepted_hits.bam

ENDUSAGE

my $version = 2.0;                    # braker.pl version number
my $alternatives_from_evidence = "true"; # output alternative transcripts based on explicit evidence from hints
my $augpath;                          # path to augustus
my $augustus_cfg_path;                # augustus config path, higher priority than $AUGUSTUS_CONFIG_PATH on system
my $augustus_bin_path;                # path to augustus folder binaries folder
my $augustus_scripts_path;	      # path to augustus scripts folder
my $AUGUSTUS_CONFIG_PATH;
if(-e $ENV{'AUGUSTUS_CONFIG_PATH'}){
  $AUGUSTUS_CONFIG_PATH = $ENV{'AUGUSTUS_CONFIG_PATH'};
}
my $AUGUSTUS_BIN_PATH;
my $AUGUSTUS_SCRIPTS_PATH;
my @bam;                              # bam file names
my $bamtools_path;                    # path to bamtools executable, higher priority than $BAMTOOLS_BIN_PATH on system
my $BAMTOOLS_BIN_PATH = $ENV{'BAMTOOLS_PATH'}; # bamtools environment variable
my @bonus;                            # array of bonus values for extrinsic file
my $bool_species = "true";            # false, if $species contains forbidden words (e.g. chmod)
my $cmdString;                        # to store shell commands
my $CPU = 1;                          # number of CPUs that can be used
my $currentDir = cwd();               # working superdirectory where program is called from
my $errorfile;                        # stores current error file name
my $errorfilesDir;                    # directory for error files
my $extrinsicCfgFile;                 # assigned extrinsic file
my $extrinsicfilesDir;                # directory for test extrinsic files
my @files;                            # contains all files in $rootDir
my $flanking_DNA;                     # length of flanking DNA, default value is min{ave. gene length/2, 10000}
my @forbidden_words;                  # words/commands that are not allowed in species name (e.g. unlink)
my $fungus = 0;                       # option for GeneMark-ET
my $gb_good_size;                     # number of LOCUS entries in 'genbank.good.gb'                         
my $genbank;                          # genbank file name
my $genemarkDir;                      # directory for GeneMark-ET output
my $GENEMARK_PATH = $ENV{'GENEMARK_PATH'}; # path to 'gmes_petap.pl' script on system, environment variable
my $GMET_path;                        # GeneMark-ET path, higher priority than $GENEMARK_PATH
my $genome;                           # name of sequence file
my $genome_length = 0;                # length of genome file
my $gff3 = 0;                         # create output file in GFF3 format
my $help;                             # print usage
my @hints;                            # input hints file names
my $hintsfile;                        # hints file later used by AUGUSTUS (RNA-Seq)
my $prot_hintsfile;                   # hints file later used by AUGUStUS (proteins)
my $limit = 10000000;                 # maximum for generic species Sp_$limit
my $logfile;                          # contains used shell commands
my @malus;                            # array of malus values for extrinsic file
my $optCfgFile;                       # optinonal extrinsic config file for AUGUSTUS
my $otherfilesDir;                    # directory for other files besides GeneMark-ET output and parameter files
my $overwrite = 0;                    # overwrite existing files (except for species parameter files)
my $parameterDir;                     # directory of parameter files for species
my $perlCmdString;                    # stores perl commands
my $printVersion = 0;                 # print version number, if set
my $SAMTOOLS_PATH = $ENV{'SAMTOOLS_PATH'}; # samtools environment variable
my $SAMTOOLS_PATH_OP;                 # path to samtools executable, higher priority than $SAMTOOLS_PATH on system
my $scriptPath=dirname($0);           # path of directory where this script is located
my $skipGeneMarkET = 0;               # skip GeneMark-ET and use provided GeneMark-ET output (e.g. from a different source)
my $skipoptimize = 0;                 # skip optimize parameter step
my $species;                          # species name
my $soft_mask = 0;                    # soft-masked flag
my $standard = 0;                     # index for standard malus/ bonus value (currently 0.1 and 1e1)optimize_augustus.pl
my $stdoutfile;                       # stores current standard output
my $string;                           # string for storing script path
my $augustus_args;                    # string that stores command line arguments to be passed to augustus
my $testsize;                         # AUGUSTUS training parameter: number of genes in a file that is
                                      # used to measure accuracy during parameter estimation with
                                      # optimize_augustus.pl. Default: 1000. If there are less than 1000
                                      # genes, all genes are used to measure accuracy. Decreasing this
                                      # parameter speeds up braker.pl but decreases prediction accuracy.
                                      # At least 300 genes are required for training AUGUSTUS.
my $useexisting = 0;                  # use existing species config and parameter files, no changes on those files
my $UTR = "off";                      # UTR prediction on/off. currently not available für new species 
my $workDir;                          # in the working directory results and temporary files are stored
my $filterOutShort;		      # filterOutShort option (see help)
my @allowedHints = ("intron");        # Currently, BRAKER only supports intron hints. Hint type from input hintsfile will
                                      # be checked.
my $crf;                              # flag that determines whether CRF training should be tried
my $nice;                             # flag that determines whether system calls should be executed with bash nice 
                                      #(default nice value)
my ($target_1, $target_2, $target_3) = 0; # variables that store AUGUSTUS accuracy after different training steps
my $prg;                              # variable to store protein alignment tool
my @prot_seq_files;                    # variable to store protein sequence file name
my @prot_aln_files;                    # variable to store protein alingment file name
my @prot_hints_files;                  # variable to store protein hints file name
my $ALIGNMENT_TOOL_PATH;              # stores path to binary of gth, spaln or exonerate for running protein alignments

# list of forbidden words for species name
@forbidden_words = ("system", "exec", "passthru", "run", "fork", "qx", "backticks", "chmod", "chown", "chroot", "unlink", "do", "eval", "kill", "rm", "mv", "grep", "cd", "top", "cp", "for", "done", "passwd", "while", "nice");  

# lists for extrinsic files
@bonus = ("1e1", "1e0", "1e2", "1e3", "1e4", "1e5");
@malus = ("0.1", "0.2", "0.4", "1.0"); 


if(@ARGV==0){
  print "$usage\n"; 
  exit(0);
}

GetOptions( 'alternatives-from-evidence=s'  => \$alternatives_from_evidence,
            'AUGUSTUS_CONFIG_PATH=s'        => \$augustus_cfg_path,
            'AUGUSTUS_BIN_PATH=s'           => \$augustus_bin_path,
            'AUGUSTUS_SCRIPTS_PATH=s'       => \$augustus_scripts_path,
	    'ALIGNMENT_TOOL_PATH=s'         => \$ALIGNMENT_TOOL_PATH,
            'bam=s'                         => \@bam,
            'BAMTOOLS_PATH=s'               => \$bamtools_path,
            'cores=i'                       => \$CPU,
            'fungus!'                       => \$fungus,
            'extrinsicCfgFile=s'            => \$extrinsicCfgFile,
            'GENEMARK_PATH=s'               => \$GMET_path,
            'genome=s'                      => \$genome,
            'gff3'                          => \$gff3,
            'rnaseq-hints=s'                => \@hints,
            'optCfgFile=s'                  => \$optCfgFile,
            'overwrite!'                    => \$overwrite,
            'SAMTOOLS_PATH=s'               => \$SAMTOOLS_PATH_OP,
            'skipGeneMark-ET!'              => \$skipGeneMarkET,
            'skipOptimize!'                 => \$skipoptimize,
            'species=s'                     => \$species,
            'softmasking!'                  => \$soft_mask,
            'testsize=i'                    => \$testsize,
            'useexisting!'                  => \$useexisting,
            'UTR=s'                         => \$UTR,
            'workingdir=s'                  => \$workDir,
	    'filterOutShort!'		    => \$filterOutShort,
	    'crf!'                          => \$crf,
	    'nice!'                         => \$nice,
            'help!'                         => \$help,
	    'prg=s'                         => \$prg,
	    'prot_seq=s'                    => \@prot_seq_files,
	    'prot_aln=s'                    => \@prot_aln_files,
	    'prot_hints=s'                  => \@prot_hints_files,
	    'augustus_args=s'               => \$augustus_args,
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
  my $tmp_dir_name = abs_path($workDir);
  $workDir = $tmp_dir_name;
}


# check the write permission of $workDir before building of the work directory
if(! -w $workDir){
  print STDERR "ERROR: Do not have write permission for $workDir.\nPlease use command 'chmod' to reset permission or specify another working directory\n"; # see autoAug.pl
  exit(1);
}

# set path to augustus config folder
if(defined($augustus_cfg_path)){
  my $last_char = substr($augustus_cfg_path, -1);
  if($last_char eq "\/"){
    chop($augustus_cfg_path);
  }
  $AUGUSTUS_CONFIG_PATH = $augustus_cfg_path;
}

# check early on whether AUGUSTUS_CONFIG_PATH is writable
if(not(-w "$AUGUSTUS_CONFIG_PATH/species")){
    print STDERR "ERROR: AUGUSTUS_CONFIG_PATH/species (in this case $AUGUSTUS_CONFIG_PATH/species) is not writeable.\n";
    exit(1)
}

# set path to augustus folders
if(defined($augustus_bin_path)){
  my $last_char = substr($augustus_bin_path, -1);
  if($last_char eq "\/"){
    chop($augustus_bin_path);
  }
  $AUGUSTUS_BIN_PATH = $augustus_bin_path;
}else{
  $AUGUSTUS_BIN_PATH = "$AUGUSTUS_CONFIG_PATH/../bin";
}

if(defined($augustus_scripts_path)){
  my $last_char = substr($augustus_scripts_path, -1);
  if($last_char eq "\/"){
    chop($augustus_scripts_path);
  }
  $AUGUSTUS_SCRIPTS_PATH = $augustus_scripts_path;
}else{
  $AUGUSTUS_SCRIPTS_PATH = "$AUGUSTUS_BIN_PATH/../scripts";
}


# set path to GeneMark-ETs gmes_petap.pl script
if(defined($GMET_path)){
  my $last_char = substr($GMET_path, -1);
  if($last_char eq "\/"){
    chop($GMET_path);
  }
  $GENEMARK_PATH = $GMET_path;
}

# set path to bamtools
if(defined($bamtools_path)){
  my $last_char = substr($bamtools_path, -1);
  if($last_char eq "\/"){
    chop($bamtools_path);
  }
  $BAMTOOLS_BIN_PATH = $bamtools_path;
}

# set path to samtools (optional)
if(defined($SAMTOOLS_PATH_OP)){
  my $last_char = substr($SAMTOOLS_PATH_OP, -1);
  if($last_char eq "\/"){
    chop($SAMTOOLS_PATH_OP);
  }
  $SAMTOOLS_PATH = $SAMTOOLS_PATH_OP;
}

# set path to alignment tool (gth, spaln or exonerate, optional)
if(defined($ALIGNMENT_TOOL_PATH)){
    my $last_char = substr($ALIGNMENT_TOOL_PATH, -1);
    if($last_char eq "\/"){
	chop($ALIGNMENT_TOOL_PATH);
    }
}

# check upfront whether any common problems will occur later # see autoAug.pl
print STDOUT "NEXT STEP: check files and settings\n"; 
check_upfront();
check_options();


# check whether RNAseq files are specified
if(!@bam && !@hints){
  print STDERR "ERROR: No RNAseq or hints files specified. Please set at least one RNAseq BAM file.\n$usage";
  exit(1);
}


# check whether bam files exist
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
    $species =~ s/\s/\_/g;
  }

  foreach my $word (@forbidden_words){
    if($species eq $word){
      print STDOUT "WARNING: $species is not allowed as a species name. ";
      $bool_species = "false";
    }
  }
}

# use standard name when no name is assigned or when it contains invalid parts
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
  if($no > $limit){
    print STDERR "ERROR: There are already $limit species folders under $AUGUSTUS_CONFIG_PATH/species/ of type 'Sp_$limit'. Please delete or move some of those folders or assign a valid species identifier with --species=name.\n";
    exit(1);
  }
  if($bool_species eq "false"){
    print STDOUT "Program will use $species instead.\n";
  }else{
    print STDOUT "No species was set. Program will use $species.\n";
  }
}
    

# check species directory
if(-d "$AUGUSTUS_CONFIG_PATH/species/$species" && !$useexisting){
  print STDERR "ERROR:$AUGUSTUS_CONFIG_PATH/species/$species already exists. Choose another species name, delete this directory or use the existing species with the option --useexisting.\n";
  exit(1);
}

if(! -d "$AUGUSTUS_CONFIG_PATH/species/$species" && $useexisting){
  print STDOUT "WARNING: $AUGUSTUS_CONFIG_PATH/species/$species does not exist. Braker will create the necessary files for species $species.\n";
  $useexisting = 0;
}

# check whether $rootDir already exists
my $rootDir = "$workDir/braker";
if (-d "$rootDir/$species" && !$overwrite){
  print STDOUT "WARNING: $rootDir/$species already exists. Braker will use existing files, if they are newer than the input files. You can choose another working directory with --workingdir=dir or overwrite it with --overwrite\n";
}

# set path and check whether assigned extrinsic file exists
if(defined($extrinsicCfgFile)){
  $extrinsicCfgFile = rel2abs($extrinsicCfgFile);
}
if(defined($extrinsicCfgFile) && ! -f $extrinsicCfgFile){
  print STDOUT "WARNING: Assigned extrinsic file $extrinsicCfgFile does not exist. Program will create extrinsic file instead.\n";
  $extrinsicCfgFile = undef;
}

# check whether genome file is set
if(!defined($genome)){
  print STDERR "ERROR: No genome file was specified. Please set a genome file.\n";
  exit(1);
}

# check whether protein sequence file is given
if(@prot_seq_files){
    @prot_seq_files = split(/[\s,]/, join(',',@prot_seq_files));
    for(my $i=0; $i<scalar(@prot_seq_files); $i++){
	if(! -f $prot_seq_files[$i]){
	    print STDERR "ERROR: protein sequence file $prot_seq_files[$i] does not exist.\n";
	    exit(1);
	}
	$prot_seq_files[$i] = rel2abs($prot_seq_files[$i]);
    }
    if(!defined($prg)){
        # if no alignment tools was specified, set Genome Threader as default
	print LOG "WARNING: No alignment tool was specified for aligning protein sequences against genome. Setting GenomeThreader as default alignment tool.\n";
	$prg="gth";
    }
}

# check whether protein alignment file is given
if(@prot_aln_files){
    @prot_aln_files = split(/[\s,]/, join(',',@prot_aln_files));
    for(my $i=0; $i<scalar(@prot_aln_files); $i++){
	if(! -f $prot_aln_files[$i]){
	    print STDERR "ERROR: protein alignment file $prot_aln_files[$i] does not exist.\n";
	    exit(1);
	}
	$prot_aln_files[$i] = rel2abs($prot_aln_files[$i]);
    }
    if(!defined($prg)){
	print STDERR "ERROR: if protein alignment file is specified, you must specify the source tool that was used to create that alignment file, i.e. --prg=gth for GenomeThreader, or --prg=spaln for Spaln2 or --prg=exonerate for Exonerate.\n";
	exit(1);
    }
}

# check whether a protein hints file is given
if(@prot_hints_files){
    @prot_hints_files = split(/[\s,]/, join(',',@prot_hints_files));
    for(my $i=0; $i<scalar(@prot_hints_files); $i++){	
	if(! -f $prot_hints_files[$i]){
	    print STDERR "ERROR: protein hints file $prot_hints_files[$i] does not exist.\n";
	    exit(1);
	}
	$prot_hints_files[$i] = rel2abs($prot_hints_files[$i]);
    }
}

# check whether alignment program is given
if(defined($prg)){
    if(not($prg=~m/gth/) and not($prg=~m/exonerate/) and not($prg=~m/spaln/)){
	print STDERR "ERROR: An alignment tool other than gth, exonerate, and spaln has been specified with option --prg=$prg. BRAKER currently only supports the options gth, exonerate and spaln.\n";
	exit(1);
    }
    if(!@prot_seq_files and !@prot_aln_files){
	print STDERR "ERROR: a protein alignment tool ($prg) has been given, but neither a protein sequence file, nor a protein alignment file generated by such a tool have been specified.\n";
	exit(1);
    }
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
  $extrinsicfilesDir = "$rootDir/$species/species/extrinsic";

  $logfile = "$otherfilesDir/braker.log";
  # create other directories if necessary
  my $bool_otherDir = "false";
  if(! -d $otherfilesDir){
    make_path($otherfilesDir);
    $bool_otherDir = "true";
  }
  open (LOG, ">>".$logfile) or die "Cannot open file: $logfile\n";
  if($bool_rootDir eq "true"){
    print LOG "\# ".(localtime).": create working directory $rootDir\n";
    print LOG "mkdir $rootDir\n\n";
  }
  if($bool_otherDir eq "true"){
    print LOG "\# ".(localtime).": create working directory $otherfilesDir\n";
    print LOG "mkdir $otherfilesDir\n\n";
  }
  if(! -d $genemarkDir){
    make_path($genemarkDir);
    print LOG "\# ".(localtime).": create working directory $genemarkDir\n";
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
    print LOG "\# ".(localtime).": create working directory $parameterDir\n";
    print LOG "mkdir $parameterDir\n\n";
  }

  if(! -d $extrinsicfilesDir){
    make_path($extrinsicfilesDir);
    print LOG "\# ".(localtime).": create working directory $extrinsicfilesDir\n";
    print LOG "mkdir $extrinsicfilesDir\n\n";
  }

  if(! -d $errorfilesDir){
    make_path($errorfilesDir);
    print LOG "\# ".(localtime).": create working directory $errorfilesDir\n";
    print LOG "mkdir $errorfilesDir\n\n";
  }

  $genome = rel2abs($genome);
  $cmdString = "cd $rootDir;";
  print LOG "\# ".(localtime).": change to working directory $rootDir\n";
  print LOG "$cmdString\n\n";
  chdir $rootDir or die ("Could not change to directory $rootDir.\n");

  new_species();                # create new species parameter files; we do this FIRST, before anything else, because if you start several processes in parallel, you might otherwise end up with those processes using the same species directory!
  check_fasta_headers($genome); # check fasta headers
  if(@prot_seq_files){
      foreach(@prot_seq_files){
	  check_fasta_headers($_);
      }
  }
###  make_rna_seq_hints();         # make hints from RNA-Seq or gtf input
  if(@prot_seq_files or @prot_aln_files or @prot_hints_files){
      make_prot_hints();
  }
###  GeneMark_ET();                # run GeneMark-ET
###  training();                   # train species-specific parameters with optimize_augustus.pl and etraining

=head
  # no extrinsic file is defined, extrinsic step is skipped or no file defined and softmasking option is used
  if(!defined($extrinsicCfgFile) || (!defined($extrinsicCfgFile) && $soft_mask)){
    extrinsic(); # use default extrinsic file
  }
  # copy extrinsic file to species parameters directory
  if(defined($extrinsicCfgFile) || (!defined($extrinsicCfgFile) && -e "$AUGUSTUS_CONFIG_PATH/species/$species/extrinsic.$species.cfg") || (!defined($extrinsicCfgFile) && -e "$AUGUSTUS_BIN_PATH/../species/$species/extrinsic.$species.cfg")){
    # using cfg files from AUGUSTUS_CONFIG_PATH has higher priority than from AUGUSTUS_BIN_PATH
    if(!defined($extrinsicCfgFile) && -e "$AUGUSTUS_CONFIG_PATH/species/$species/extrinsic.$species.cfg"){
      $extrinsicCfgFile = "$AUGUSTUS_CONFIG_PATH/species/$species/extrinsic.$species.cfg";
    }elsif(!defined($extrinsicCfgFile) && -e "$AUGUSTUS_BIN_PATH/../species/$species/extrinsic.$species.cfg"){
      $extrinsicCfgFile = "$AUGUSTUS_BIN_PATH/../species/$species/extrinsic.$species.cfg";
    }
    @_ = split(/\//, $extrinsicCfgFile);
    if(!uptodate([$extrinsicCfgFile],["$parameterDir/$species/".$_[-1]])  || $overwrite){
	print STDOUT "NEXT STEP: copy extrinsic file to working directory\n";
	$cmdString = "";
	if($nice){
	    $cmdString .= "nice ";
	}
	$cmdString .= "cp $extrinsicCfgFile $parameterDir/$species/$_[-1]";
	print LOG "\# ".(localtime).": copy extrinsic file to working directory\n";
	print LOG "$cmdString\n\n";
	system("$cmdString")==0 or die("failed to execute: $cmdString!\n");
	print STDOUT "extrinsic file copied.\n";
    }
  }

###  augustus(); # run augustus
  if(!uptodate(["$otherfilesDir/augustus.gff"],["$otherfilesDir/augustus.aa"])  || $overwrite){
###    getAnnoFasta("$otherfilesDir/augustus.gff"); # create protein sequence file
  }
  if(!uptodate(["$otherfilesDir/augustus.gff"],["$otherfilesDir/augustus.gtf"])  || $overwrite){
###    make_gtf("$otherfilesDir/augustus.gff"); # convert output to gtf and gff3 (if desired) format
  }
###  clean_up(); # delete all empty files
=cut
  close(LOG) or die("Could not close log file $logfile!\n");
}



                           ############### sub functions ##############


         ####################### make hints #########################
# make hints from BAM files and optionally combine it with additional hints file
sub make_rna_seq_hints{
  my $bam_hints;
  my $hintsfile_temp = "$otherfilesDir/hintsfile.temp.gff";
  $hintsfile = "$otherfilesDir/hintsfile.gff";
  # from RNA-Seq data in bam format
  if(@bam){
    my $bam_temp = "$otherfilesDir/bam2hints.temp.gff";
    for(my $i=0; $i<scalar(@bam); $i++){
      $errorfile = "$errorfilesDir/bam2hints.$i.stderr";
      if(!uptodate([$bam[$i]],[$hintsfile])  || $overwrite){
        $bam[$i] = check_bam_headers($bam[$i]);
        print STDOUT "NEXT STEP: make hints from BAM file $bam[$i]\n";
        if(-e "$AUGUSTUS_CONFIG_PATH/../bin/bam2hints"){
          $augpath = "$AUGUSTUS_CONFIG_PATH/../bin/bam2hints";
        }else{
          $augpath = "$AUGUSTUS_BIN_PATH/bin/bam2hints";
        }
	$cmdString = "";
	if($nice){
	    $cmdString .= "nice ";
	}
	$cmdString .= "$augpath --intronsonly --in=$bam[$i] --out=$bam_temp 2>$errorfile";
        print LOG "\# ".(localtime).": make hints from BAM file $bam[$i]\n";
        print LOG "$cmdString\n\n";
        system("$cmdString")==0 or die("failed to execute: $cmdString!\n");
	$cmdString = "";
	if($nice){
	    $cmdString .= "nice ";
	}
        $cmdString .= "cat $bam_temp >>$hintsfile_temp";
        print LOG "\# ".(localtime).": add hints from BAM file $bam[$i] to hints file\n";
        print LOG "$cmdString\n\n";
        system("$cmdString")==0 or die("failed to execute: $cmdString!\n");
        print STDOUT "hints from BAM file $bam[$i] added.\n";
      }
    }
    unlink($bam_temp);
  }
  # from gtf input files
  if(@hints){
    for(my $i=0; $i<scalar(@hints); $i++){
      if(!uptodate([$hints[$i]],[$hintsfile]) || $overwrite){
	  pricnt STDOUT "NEXT STEP: add hints from file $hints[$i]\n";
	  $cmdString = "";
	  if($nice){
	      $cmdString .= "nice ";
	  }
	  $cmdString .= "cat $hints[$i] >> $hintsfile_temp";
	  print LOG "\# ".(localtime).": add hints from file $hints[$i]\n";
	  print LOG "$cmdString\n\n";
	  system("$cmdString")==0 or die("failed to execute: $cmdString!\n");
      }
    }
  }
  if(-f $hintsfile_temp || $overwrite){
    if(!uptodate([$hintsfile_temp],[$hintsfile]) || $overwrite){
	join_mult_hints($hintsfile_temp, "rnaseq");
    }
    if(!uptodate([$hintsfile_temp],[$hintsfile]) || $overwrite){
      print STDOUT "NEXT STEP: filter introns, find strand and change score to \'mult\' entry\n";
      $string = find("filterIntronsFindStrand.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
      $errorfile = "$errorfilesDir/filterIntronsFindStrand.stderr";
      $perlCmdString = "";
      if($nice){
	  $perlCmdString .= "nice ";
      }
      $perlCmdString .= "perl $string $genome $hintsfile_temp --score 1>$hintsfile 2>$errorfile";
      print LOG "\# ".(localtime).": filter introns, find strand and change score to \'mult\' entry\n";
      print LOG "$perlCmdString\n\n";
      system("$perlCmdString")==0 or die("failed to execute: $perlCmdString!\n");
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

sub make_prot_hints{
    my $prot_hints;
    my $prot_hints_file_temp = "$otherfilesDir/prot_hintsfile.temp.gff";
    $prot_hintsfile = "$otherfilesDir/prot_hintsfile.gff";
    my $alignment_outfile = "$otherfilesDir/protein_alignment_$prg.gff3";   
    # change to working directory
    $cmdString = "cd $otherfilesDir";
    print LOG "\# ".(localtime).": Changing to $otherfilesDir\n";
    print LOG "$cmdString\n";
    chdir $otherfilesDir or die ("Failed to execude $cmdString!\n");
    # from fasta files
    if(@prot_seq_files){
	print STDOUT "NEXT STEP: running alignment tool $prg\n";
	$string = find("startAlign.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
	$errorfile = "$errorfilesDir/startAlign.stderr";
	$logfile = "$errorfilesDir/startAlign.stdout";
	for(my $i=0; $i<scalar(@prot_seq_files); $i++){
	    if(!uptodate([$prot_seq_files[$i]],[$prot_hintsfile])  || $overwrite){
		$perlCmdString = "";
		if($nice){
		    $perlCmdString .= "nice ";
		}
		$perlCmdString .= "perl $string --genome=$genome --prot=$prot_seq_files[$i] --ALIGNMENT_TOOL_PATH=$ALIGNMENT_TOOL_PATH ";
		if($prg eq "gth"){
		    $perlCmdString .= "--prg=gth ";
		    print LOG "\# ".(localtime).": running Genome Threader to produce protein to genome alignments\n";
		}elsif($prg eq "exonerate"){
		    $perlCmdString .= "--prg=exonerate ";
		    print LOG "\# ".(localtime).": running Exonerate to produce protein to genome alignments\n";
		}elsif($prg eq "spaln"){
		    $perlCmdString .= "--prg=spaln ";
		    print LOG "\# ".(localtime).": running Spaln to produce protein to genome alignments\n";
		}
		if($CPU>1){
		    $perlCmdString .= "--CPU=$CPU ";
		}
		if($nice){
		    $perlCmdString .= "--nice ";
		}
		$perlCmdString .= ">> $logfile 2>>$errorfile";
		print LOG "$perlCmdString\n\n";
		system("$perlCmdString")==0 or die ("failed to execute: $perlCmdString!\n");
		print STDOUT "Alignments from file $prot_seq_files[$i] created.\n";
		$cmdString = "cat align_$prg/$prg.concat.aln >> $alignment_outfile";	    
		print LOG "\# ".(localtime).": concatenating alignment file to $alignment_outfile\n";
		print LOG "$cmdString\n\n";
		system("$cmdString")==0 or die("Failed to execute $cmdString!\n");
		print LOG "\# ".(localtime).": moving startAlign output files\n";
		$cmdString = "mv startAlign_$prg.log startAlign_$prg.log$i";
		print LOG "$cmdString\n";
		system("$cmdString")==0 or die("Failed to execute $cmdString!\n");
		$cmdString = "mv tmp_$prg tmp_$prg$i";
		print LOG "$cmdString\n";
		system("$cmdString")==0 or die("Failed to execute $cmdString!\n");
		$cmdString = "mv align_$prg align_$prg$i";
		print LOG "$cmdString\n";
		system("$cmdString")==0 or die("Failed to execute $cmdString!\n\n");	  
	    }else{
		print STDOUT "Skipping running alignment tool because files $prot_seq_files[$i] and $prot_hintsfile were up to date.\n";		
	    }  
	}
    }
    # convert pipeline created protein alignments to protein hints
    if(@prot_seq_files && -e $alignment_outfile){
	if(!uptodate([$alignment_outfile], [$prot_hintsfile]) || $overwrite){
	    aln2hints($alignment_outfile, $prot_hints_file_temp);
	}
    }
    # convert command line specified protein alignments to protein hints
    if(@prot_aln_files){
	for(my $i=0; $i<scalar(@prot_aln_files); $i++){
	    if(!uptodate([$prot_aln_files[$i]], [$prot_hintsfile]) || $overwrite){
		aln2hints($prot_aln_files[$i], $prot_hints_file_temp);
	    }else{
		print "Skipped converting alignment file $prot_aln_files[$i] to hints because it was up to date with $prot_hintsfile\n";
	    }
	}
    }
    # append command line specified hints
    if(@prot_hints_files){
	for(my $i=0; $i<scalar(@prot_hints_files); $i++){
	    if(!uptodate([$prot_hints_files[$i]], [$prot_hintsfile]) || $overwrite){
		print STDOUT "NEXT STEP: appending hints file $prot_hints_files[$i] to protein hints file $prot_hints_file_temp\n";
		print LOG "\# ".(localtime).": Appending hints file $prot_hints_files[$i] to protein hints file $prot_hints_file_temp\n";
		$cmdString = "cat $prot_hints_files[$i] >> $prot_hints_file_temp";
		print LOG $cmdString."\n";
		system($cmdString)==0 or die ("Could not execute $cmdString!\n");
		print STDOUT "File appended.\n";
	    }
	}
    }
    if(-f $prot_hints_file_temp || $overwrite){
	if(!uptodate([$prot_hints_file_temp],[$prot_hintsfile])|| $overwrite){
	    join_mult_hints($prot_hints_file_temp, "prot");
	    print STDOUT "NEXT STEP: moving $prot_hints_file_temp to $prot_hintsfile\n";
	    print LOG "\# ".(localtime).": moving $prot_hints_file_temp to $prot_hintsfile\n";
$cmdString = "mv $prot_hints_file_temp $prot_hintsfile";
	    print LOG "$cmdString\n";
	    system($cmdString)==0 or die ("Could not execute $cmdString!\n");
	    print LOG "Deleting $prot_hints_file_temp\n";
	    unlink($prot_hints_file_temp);
	    print STDOUT "file moved.\n";
	}
    }
    if(-z $prot_hintsfile){
	print STDERR "ERROR: The hints file is empty. There were no protein alignments.\n";
	exit(1);
    }
}


sub aln2hints{
    my $aln_file = shift;
    if(! (-z $aln_file)){
	my $out_file_name = "$otherfilesDir/prot_hintsfile.aln2hints.temp.gff";
	my $final_out_file = shift;
	print STDOUT "NEXT STEP: converting alignment file $aln_file to protein hints file\n";
	print LOG "\# ".(localtime).": Converting protein alignment file $aln_file to hints for AUGUSTUS\n";
	$perlCmdString = "perl ";
	if($nice){
	    $perlCmdString .= "nice ";
	}
	$string = find("align2hints.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
	$perlCmdString .= "$string --in=$aln_file --out=$out_file_name ";
	if($prg eq "spaln"){
	    $perlCmdString .= "--prg=spaln";
	}elsif($prg eq "gth"){
	    $perlCmdString .= "--prg=gth";
	}elsif($prg eq "exonerate"){
	    $perlCmdString .= "--prg=exonerate --genome_file=$genome";
	}
	print LOG "$perlCmdString\n";
	system("$perlCmdString")==0 or die ("Failed to execute $perlCmdString!\n\n");
	print STDOUT "$out_file_name protein hints file created.\n";
	$cmdString = "cat $out_file_name >> $final_out_file";
	print LOG "\# ".(localtime).": concatenating protein hints from $out_file_name to $final_out_file\n";
	print LOG $cmdString."\n";
	system("$cmdString")==0 or die ("Could not execute $cmdString!\n");
    }else{
	print "Alignment file $aln_file was empty!\n";
	print LOG "\# ".(localtime).": Alignment file $aln_file was empty!\n";
    }
}

sub join_mult_hints{
    my $hints_file_temp = shift;
    my $type = shift; # rnaseq or prot or whatever will follow
    my $hintsfile_temp_sort = "$otherfilesDir/hints.$type.temp.sort.gff";
    print STDOUT "NEXT STEP: sort hints\n";
    $cmdString = "";
    if($nice){
	$cmdString .= "nice ";
    }
    $cmdString .= "cat $hints_file_temp | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -n -k 3,3 | sort -s -k 1,1 >$hintsfile_temp_sort";
    print LOG "\# ".(localtime).": sort hints of type $type\n";
    print LOG "$cmdString\n\n";
    system("$cmdString")==0 or die("failed to execute: $cmdString!\n");
    print STDOUT "hints sorted.\n";
    print STDOUT "NEXT STEP: summarize multiple identical hints of type to one\n";
    $string = find("join_mult_hints.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
    $errorfile = "$errorfilesDir/join_mult_hints.$type.stderr";
    $perlCmdString = "";
    if($nice){
	$perlCmdString .= "nice ";
    }
    $perlCmdString .= "perl $string <$hintsfile_temp_sort >$hints_file_temp 2>$errorfile";
    print LOG "\# ".(localtime).": join multiple hints\n";
    print LOG "$perlCmdString\n\n";
    system("$perlCmdString")==0 or die("failed to execute: $perlCmdString!\n");
    print STDOUT "hints joined.\n";
    unlink($hintsfile_temp_sort);
}
    



         ####################### GeneMark-ET #########################
# start GeneMark-ET and convert its output to real gtf format
sub GeneMark_ET{
  if(!$skipGeneMarkET){
    if(!uptodate([$genome,$hintsfile],["$genemarkDir/genemark.gtf"])  || $overwrite){
      $cmdString = "cd $genemarkDir";
      print LOG "\# ".(localtime).": change to GeneMark-ET directory $genemarkDir\n";
      print LOG "$cmdString\n\n";
      chdir $genemarkDir or die ("Could not change to directory $genemarkDir.\n");

      print STDOUT "NEXT STEP: execute GeneMark-ET\n"; 
      $string = "$GENEMARK_PATH/gmes_petap.pl";
      $errorfile = "$errorfilesDir/GeneMark-ET.stderr";
      $stdoutfile = "$otherfilesDir/GeneMark-ET.stdout";
      $perlCmdString = "";
      if($nice){
	  $perlCmdString .= "nice ";
      }
      $perlCmdString .= "perl $string --sequence=$genome --ET=$hintsfile --cores=$CPU";
      if($fungus){
        $perlCmdString .= " --fungus";
      }
      if($soft_mask){
      #  $perlCmdString .= " --soft_mask"; # version prior to 4.29
        $perlCmdString .= " --soft 1000"; # version 4.29
      }
      $perlCmdString .= " 1>$stdoutfile 2>$errorfile";
      print LOG "\# ".(localtime).": execute GeneMark-ET\n";
      print LOG "$perlCmdString\n\n";
      system("$perlCmdString")==0 or die("failed to execute: $perlCmdString\n");
      print STDOUT "GeneMark-ET finished.\n";
      $cmdString = "";
      if($nice){
	  $cmdString .= "nice ";
      }
      $cmdString .= "cd $rootDir";
      print LOG "\# ".(localtime).": change to working directory $rootDir\n";
      print LOG "$cmdString\n\n";
      chdir $rootDir or die ("Could not change to directory $rootDir.\n");
    }
  }

  # convert GeneMark-ET output to gtf format with doublequotes (for older GeneMark-ET versions) and filter genes for training
  if(!uptodate(["$genemarkDir/genemark.gtf", $hintsfile],["$genemarkDir/genemark.c.gtf","$genemarkDir/genemark.f.good.gtf", "$genemarkDir/genemark.average_gene_length.out"])  || $overwrite){
    print STDOUT "NEXT STEP: convert GeneMark-ET to real gtf format\n"; 
    $string=find("filterGenemark.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
    $errorfile = "$errorfilesDir/filterGenemark.stderr";
    $stdoutfile = "$otherfilesDir/filterGenemark.stdout";
    $perlCmdString = "";
    if($nice){
	$perlCmdString .= "nice ";
    }
    if(!$filterOutShort){
        $perlCmdString .= "perl $string --genemark=$genemarkDir/genemark.gtf --introns=$hintsfile 1>$stdoutfile 2>$errorfile";
    }else{
	print STDOUT "Option activated: Filtering out training genes from GeneMark that are too short (upstream intron)\n";
	$perlCmdString .= "perl $string --genemark=$genemarkDir/genemark.gtf --introns=$hintsfile --filterOutShort 1>$stdoutfile 2>$errorfile";
    }
    print LOG "\# ".(localtime).": convert GeneMark-ET output to real gtf format\n";
    print LOG "$perlCmdString\n\n";
    system("$perlCmdString")==0 or die("failed to execute: $perlCmdString!\n");
    print STDOUT "GeneMark-ET conversion.\n";
  }
}



         ####################### create a new species #########################
# create a new species $species and extrinsic file from generic one 
sub new_species{
  $augpath = "$AUGUSTUS_CONFIG_PATH/species/$species";
  if((!uptodate([$augpath."/$species\_metapars.cfg"],[$augpath."/$species\_parameters.cfg", $augpath."/$species\_exon_probs.pbl"]) && !$useexisting) || ! -d "$AUGUSTUS_CONFIG_PATH/species/$species"){
      if(-d "$AUGUSTUS_CONFIG_PATH/species"){
	  if(-w "$AUGUSTUS_CONFIG_PATH/species"){
	      print STDOUT "NEXT STEP: create new species $species in $AUGUSTUS_CONFIG_PATH/species\n";
      
	      $string=find("new_species.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
	      $errorfile = "$errorfilesDir/new_species.stderr";
	      $perlCmdString = "";
	      if($nice){
		  $perlCmdString .= "nice ";
	      }
	      $perlCmdString .= "perl $string --species=$species --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH 1> /dev/null 2>$errorfile";
	      print LOG "\# ".(localtime).": create new species $species\n";
	      print LOG "$perlCmdString\n\n";
	      system("$perlCmdString")==0 or die("Failed to create new species with new_species.pl, check write permissions in $AUGUSTUS_CONFIG_PATH/species directory! $!\n");
	  }else{
	      print STDERR "Directory $AUGUSTUS_CONFIG_PATH/species is not writable! You must make the directory AUGUSTUS_CONFIG_PATH/species writable or specify another AUGUSTUS_CONFIG_PATH!\n";
	  }
      }else{
	  print STDERR "Directory $AUGUSTUS_CONFIG_PATH/species does not exist. Please check that AUGUSTUS_CONFIG_PATH is set, correctly; and/or create species directory!\n";
      }
  }
}


# create default extrinsic file (without BLAST)

sub extrinsic{
    my $extrinsic;
  if(-e "$AUGUSTUS_CONFIG_PATH/extrinsic/extrinsic.M.RM.E.W.cfg"){
    $extrinsic = "$AUGUSTUS_CONFIG_PATH/extrinsic/extrinsic.M.RM.E.W.cfg";
    print STDOUT "Will use $AUGUSTUS_CONFIG_PATH/extrinsic/extrinsic.M.RM.E.W.cfg as template for this project's extrinsic.cfg\n";
  }elsif(-e "$AUGUSTUS_BIN_PATH/../config/extrinsic/extrinsic.M.RM.E.W.cfg"){
    $extrinsic = "$AUGUSTUS_BIN_PATH/../config/extrinsic/extrinsic.M.RM.E.W.cfg";
    print STDOUT "Will use $AUGUSTUS_BIN_PATH/../config/extrinsic/extrinsic.M.RM.E.W.cfg as template for this project's extrinsic.cfg\n";
  }else{
      die "Cannot find extrinsic template file $AUGUSTUS_CONFIG_PATH/extrinsic/extrinsic.M.RM.E.W.cfg or $AUGUSTUS_BIN_PATH/../config/extrinsic/extrinsic.M.RM.E.W.cfg. Please check whether AUGUSTUS_CONFIG_PATH (and/or AUGUSTUS_BIN_PATH) are correct. Looking for template file in extrinsic folder of AUGUSTUS_CONFIG_PATH and relative ../config/extrinsic to AUGUSTUS_BIN_PATH.\n";
  }
  my $extrinsic_cp = "$AUGUSTUS_CONFIG_PATH/species/$species/extrinsic.$species.cfg";
  if(!defined($extrinsicCfgFile) || (defined($extrinsicCfgFile) && ! -e $extrinsicCfgFile) ){
    if(!defined($extrinsicCfgFile)){
      print STDOUT "No extrinsic file assigned. Program will create one.\n"; # other check/warning see beginning 
    }
    if(! -e $extrinsic_cp){
      print STDOUT "NEXT STEP: create extrinsic file: $extrinsic_cp\n";
      print LOG "\# ".(localtime).": create extrinsic file\n";
      print STDOUT "species $species created.\n";

      open (EXTRINSIC, $extrinsic) or die "Cannot open file: $extrinsic\n";
      open (OUT, ">".$extrinsic_cp) or die "Cannot open file: $extrinsic_cp\n";
      my $GENERAL = "false";
      my $bonus_idx = $standard; # 0 => 1e1  (currently standard)
      my $malus_idx = $standard; # 0 => 0.1  (currently standard)
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
          print OUT "   exonpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print OUT "       exon        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print OUT " intronpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print OUT "     intron        1      $malus[$malus_idx]  M    1  1e+100  RM  1  1.15    E 1  $bonus[$bonus_idx]    W 1    1\n";
          print OUT "    CDSpart        1   1    1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print OUT "        CDS        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print OUT "    UTRpart        1   1    1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print OUT "        UTR        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print OUT "     irpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print OUT "nonexonpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print OUT "  genicpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";

          print LOG "\# ".(localtime).": edit extrinsic file and add\n$_\n";
          print LOG "\n";
          print LOG "      start        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print LOG "       stop        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print LOG "        tss        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print LOG "        tts        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print LOG "        ass        1   1    1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print LOG "        dss        1   1    1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print LOG "   exonpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print LOG "       exon        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print LOG " intronpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print LOG "     intron        1      $malus[$malus_idx]  M    1  1e+100  RM  1  1.15    E 1  $bonus[$bonus_idx]    W 1    1\n";
          print LOG "    CDSpart        1   1    1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print LOG "        CDS        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print LOG "    UTRpart        1   1    1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print LOG "        UTR        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print LOG "     irpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print LOG "nonexonpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
          print LOG "  genicpart        1        1  M    1  1e+100  RM  1     1    E 1    1    W 1    1\n";
        }
      }
      close(OUT) or die("Could not close extrinsic file $extrinsic_cp!\n");
      close(EXTRINSIC) or die("Could not close extrinsic file $extrinsic!\n");
      print STDOUT "extrinsic file created.\n";
      $extrinsicCfgFile = $extrinsic_cp;
    }else{
      $extrinsicCfgFile = $extrinsic_cp;
    }
  }else{
    print STDOUT "No extrinsic file was created. Program uses assigned extrinsic file: $extrinsicCfgFile\n";
  }
}

         ####################### train AUGUSTUS #########################
# create genbank file and train AUGUSTUS parameters for given species
sub training{
  # get average gene length for flanking region
  print STDOUT "NEXT STEP: get flanking DNA size\n"; 
  # compute average gene length from genemark.gtf                                                       
  # The following code was introduced after Katharina observed a negative average gene length in the GeneMark-ET output file. However, this seems to be an artefact from an interrupted GeneMark-ET run. Therefore, we return to the original parsing of GeneMark-ET output file.
  #    open (GENEMARKGTF, "<", "$genemarkDir/genemark.gtf") or die("Cannot open file: $genemarkDir/genemark.gtf\n");
  #  my $cumulativeGeneMarkGeneLength = 0;
  #  my $geneMarkGenesNo = 0;
  #  my %geneMarkGeneHash;
  #  while(<GENEMARKGTF>){
  #      if(m/\tCDS\t/){
  #	  @_ = split(/\t/);
  #	  $cumulativeGeneMarkGeneLength += $_[4]-$_[3]+1;
  #	  if(not(defined($geneMarkGeneHash{$_[8]}))){$geneMarkGeneHash{$_[8]} = 1;}
  #      }
  #  }
  #  $geneMarkGenesNo = keys %geneMarkGeneHash;
  #  my $average_length = $cumulativeGeneMarkGeneLength/$geneMarkGenesNo;
  #  close (GENEMARKGTF) or die("Cannot close file: $genemarkDir/genemark.gtf\n");
  open (GENELENGTH, "<$genemarkDir/genemark.average_gene_length.out") or die "Cannot open file: $genemarkDir/genemark.average_gene_length.out\n";
    @_ = split(/\t/,<GENELENGTH>);
    my $average_length = $_[0];
    close(GENELENGTH) or die("Could not close file $genemarkDir/genemark.average_gene_length.out!\n");
    if($average_length < 0){
        die("Average gene length of genes predicted by GeneMark-ET is negative. This indicates GeneMark-ET did not finish, properly. Please delete the GeneMark-ET folder and restart BRAKER1!\n");
    }
    $flanking_DNA = min((floor($average_length/2), 10000));
    if($flanking_DNA < 0){ # added by Katharina Hoff
        print STDOUT "\$flanking_DNA has the value $flanking_DNA , which is smaller than 0. Something must have gone wrong, there. Replacing by value 500. It is completely unclear whether 500 is a good or a bad choice.\n"; # added by Katharina Hoff
        $flanking_DNA = 500; # added by Katharina Hoff
    }
    # create genbank file from fasta input and GeneMark-ET output
    $string = find("gff2gbSmallDNA.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
    $genbank = "$otherfilesDir/genbank.gb";
    if(!uptodate([$genome,"$genemarkDir/genemark.c.gtf"],[$genbank])  || $overwrite){
        print STDOUT "NEXT STEP: create genbank file\n";
        $errorfile = "$errorfilesDir/gff2gbSmallDNA.stderr";
        if(-z "$genemarkDir/genemark.c.gtf"){
            print STDERR "ERROR: The GeneMark-ET output file is empty!\n";
            exit(1);
        }
        $perlCmdString = "";
        if($nice){
	        $perlCmdString .= "nice ";
        }
        $perlCmdString .= "perl $string $genemarkDir/genemark.c.gtf $genome $flanking_DNA $genbank 2>$errorfile";
        print LOG "\# ".(localtime).": create genbank file\n";
        print LOG "$perlCmdString\n\n";
        system("$perlCmdString")==0 or die("failed to execute: $perlCmdString!\n");
        print STDOUT "genbank file created.\n";  
    }

    # filter genbank file and only use the genes considered "good" (i.e. genes whose introns are represented in hints file) 
    if(!uptodate([$genbank, "$genemarkDir/genemark.f.good.gtf"],["$otherfilesDir/genbank.good.gb"])  || $overwrite){
        print STDOUT "NEXT STEP: filter genbank file\n";
        $string=find("filterGenesIn_mRNAname.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
        $errorfile = "$errorfilesDir/filterGenesIn_mRNAname.stderr";
        if(-z "$genemarkDir/genemark.f.good.gtf"){
            print STDERR "ERROR: The GeneMark-ET output file does not contain genes that are supported by the intron hints.\n";
            exit(1);
        }
        $perlCmdString = "";
        if($nice){
	        $perlCmdString .= "nice ";
        }
        $perlCmdString .= "perl $string $genemarkDir/genemark.f.good.gtf $genbank 1>$otherfilesDir/genbank.good.gb 2>$errorfile";
        print LOG "\# ".(localtime).": filter genbank file\n";
        print LOG "$perlCmdString\n\n";
        system("$perlCmdString")==0 or die("failed to execute: $perlCmdString!\n");
        print STDOUT "genbank file filtered.\n";
    }

    # split into training and test set
    if(!uptodate(["$otherfilesDir/genbank.good.gb"],["$otherfilesDir/genbank.good.gb.test", "$otherfilesDir/genbank.good.gb.train"])  || $overwrite){
        print STDOUT "NEXT STEP: split genbank file into train and test file\n";
        $string = find("randomSplit.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
        $errorfile = "$errorfilesDir/randomSplit.stderr";
        if($nice){
            $gb_good_size = `nice grep -c LOCUS $otherfilesDir/genbank.good.gb`;
        }else{
	        $gb_good_size = `grep -c LOCUS $otherfilesDir/genbank.good.gb`;
        }
        if($gb_good_size < 300){
            print STDOUT "WARNING: Number of good genes is low ($gb_good_size). Recomended are at least 300 genes\n";
        }
        if($gb_good_size == 0){
            print STDERR "ERROR: Number of good genes is 0, so the parameters cannot be optimized. Recomended are at least 300 genes\n";
            exit(1);
        }
        if($gb_good_size > 1000){
	        $testsize = 1000;
	        $perlCmdString = "";
	        if($nice){
	            $perlCmdString .= "nice ";
	        }
	        $perlCmdString .= "perl $string $otherfilesDir/genbank.good.gb $testsize 2>$errorfile";
	        print LOG "\# ".(localtime).": split genbank file into train and test file\n";
	        print LOG "$perlCmdString\n\n";
	        system("$perlCmdString")==0 or die("failed to execute: $perlCmdString!\n");
	        print STDOUT "genbank file splitted.\n";
        }
    }

    if(!$useexisting){
        # train AUGUSTUS for the first time
        if(!uptodate(["$otherfilesDir/genbank.good.gb.train","$otherfilesDir/genbank.good.gb"],["$otherfilesDir/firstetraining.stdout"])){
            # set "stopCodonExcludedFromCDS" to true
            print STDOUT "NEXT STEP: Setting value of \"stopCodonExcludedFromCDS\" in $AUGUSTUS_CONFIG_PATH/species/$species/$species\_parameters.cfg to \"true\"\n"; # see autoAugTrain.pl
            print LOG "\# ".(localtime).": Setting value of \"stopCodonExcludedFromCDS\" in $AUGUSTUS_CONFIG_PATH/species/$species/$species\_parameters.cfg to \"true\"\n";
            setParInConfig($AUGUSTUS_CONFIG_PATH."/species/$species/$species\_parameters.cfg", "stopCodonExcludedFromCDS", "true"); # see autoAugTrain.pl

            # first try with etraining
            print STDOUT "NEXT STEP: first etraining\n"; 
            $augpath = "$AUGUSTUS_BIN_PATH/etraining";
            $errorfile = "$errorfilesDir/firstetraining.stderr";
            $stdoutfile = "$otherfilesDir/firstetraining.stdout";
            if($nice){
                $gb_good_size = `nice grep -c LOCUS $otherfilesDir/genbank.good.gb`;
            }else{
	            $gb_good_size = `grep -c LOCUS $otherfilesDir/genbank.good.gb`;
            }
            $cmdString = "";
            if($nice){
	            $cmdString .= "nice ";
            }
            if($gb_good_size <= 1000){
                $cmdString .= "$augpath --species=$species --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH $otherfilesDir/genbank.good.gb 1>$stdoutfile 2>$errorfile";
                $testsize = $gb_good_size - 1;
            }else{
                $testsize = 1000;
                $cmdString .= "$augpath --species=$species --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH $otherfilesDir/genbank.good.gb.train 1>$stdoutfile 2>$errorfile";
            }
            print LOG "\# ".(localtime).": first etraining\n";
            print LOG "$cmdString\n\n";
            system("$cmdString")==0 or die("failed to execute first etraining: $!\n");
            print STDOUT "first training complete.\n";

            # set "stopCodonExcludedFromCDS" to false and run etraining again if necessary
            my $t_b_t = $gb_good_size - $testsize;      # see autoAugTrain.pl
            my $err_stopCodonExcludedFromCDS;
            if($nice){
	            $err_stopCodonExcludedFromCDS = `nice grep -c "exon doesn't end in stop codon" $errorfile`;
            }else{
	            $err_stopCodonExcludedFromCDS = `grep -c "exon doesn't end in stop codon" $errorfile`; # see autoAugTrain.pl
            }
            my $err_rate =  $err_stopCodonExcludedFromCDS / $t_b_t;  # see autoAugTrain.pl
            print STDOUT "Error rate of missing stop codon is $err_rate\n";  # see autoAugTrain.pl
            print LOG "\# ".(localtime).": Error rate of missing stop codon is $err_rate\n"; # see autoAugTrain.pl
            if($err_rate >= 0.5){ # see autoAugTrain.pl
                print STDOUT "The appropriate value for \"stopCodonExcludedFromCDS\" seems to be \"false\".\n"; # see autoAugTrain.pl
                print STDOUT "next step: Setting value of \"stopCodonExcludedFromCDS\" in $AUGUSTUS_CONFIG_PATH/species/$species/$species\_parameters.cfg to \"false\"\n"; # see autoAugTrain.pl
                setParInConfig($AUGUSTUS_CONFIG_PATH."/species/$species/$species\_parameters.cfg", "stopCodonExcludedFromCDS", "false");  # see autoAugTrain.pl
                print STDOUT "NEXT STEP: Trying etraining again\n";
                print LOG "\# ".(localtime).": Trying etraining again\n";
                print LOG "$cmdString\n\n";
                system("$cmdString")==0 or die("failed to execute second etraining: $!\n");
                print STDOUT "trying etraining again complete.\n";
            }

            # adjust the stop-codon frequency in species_parameters.cfg according to train.out
            print LOG "\# ".(localtime).": adjust the stop-codon frequency in species_parameters.cfg according to $stdoutfile\n";
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
            print LOG "\# ".(localtime).": Setting frequency of stop codons to tag=$freqOfTag, taa=$freqOfTaa, tga=$freqOfTga.\n";
            print STDOUT "NEXT STEP: Setting frequency of stop codons to tag=$freqOfTag, taa=$freqOfTaa, tga=$freqOfTga.\n";
            setParInConfig($AUGUSTUS_CONFIG_PATH."/species/$species/$species\_parameters.cfg", "/Constant/amberprob", $freqOfTag);  # see autoAugTrain.pl
            setParInConfig($AUGUSTUS_CONFIG_PATH."/species/$species/$species\_parameters.cfg", "/Constant/ochreprob", $freqOfTaa);  # see autoAugTrain.pl
            setParInConfig($AUGUSTUS_CONFIG_PATH."/species/$species/$species\_parameters.cfg", "/Constant/opalprob", $freqOfTga);   # see autoAugTrain.pl
            print STDOUT "frequency adjusted\n";
        }

        # first test
        if(!uptodate(["$otherfilesDir/genbank.good.gb.test", "$otherfilesDir/genbank.good.gb"],["$otherfilesDir/firsttest.stdout"])  || $overwrite){
            print STDOUT "NEXT STEP: first test\n";
            $augpath = "$AUGUSTUS_BIN_PATH/augustus";
            $errorfile = "$errorfilesDir/firsttest.stderr";
            $stdoutfile = "$otherfilesDir/firsttest.stdout";
            if($nice){
                $gb_good_size = `nice grep -c LOCUS $otherfilesDir/genbank.good.gb`;
            }else{
	            $gb_good_size = `grep -c LOCUS $otherfilesDir/genbank.good.gb`;
            }
            $cmdString = "";
            if($nice){
	            $cmdString .= "nice ";
            }
            if($gb_good_size <= 1000){
                $cmdString .= "$augpath --species=$species --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH $otherfilesDir/genbank.good.gb 1>$stdoutfile 2>$errorfile";
            }else{
                $cmdString .= "$augpath --species=$species --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH $otherfilesDir/genbank.good.gb.test 1>$stdoutfile 2>$errorfile";
            }
            print LOG "\# ".(localtime).": first test\n";
            print LOG "$cmdString\n\n";
            system("$cmdString")==0 or die("failed to execute: $cmdString!\n");
	        $target_1 = accuracy_calculator($stdoutfile);
	        print LOG "\# The accuracy after initial training (no optimize_augustus.pl, no CRF) is $target_1\n";      
            print STDOUT "first test finished.\n";
        }

        # optimize parameters
        if(!$skipoptimize){
            if(!uptodate(["$otherfilesDir/genbank.good.gb.train","$otherfilesDir/genbank.good.gb.test", "$otherfilesDir/genbank.good.gb"],[$AUGUSTUS_CONFIG_PATH."/species/$species/$species\_exon_probs.pbl", $AUGUSTUS_CONFIG_PATH."/species/$species/$species\_parameters.cfg", $AUGUSTUS_CONFIG_PATH."/species/$species/$species\_weightmatrix.txt"])){
                print STDOUT "NEXT STEP: optimize AUGUSTUS parameter\n";
                $string=find("optimize_augustus.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
                $errorfile = "$errorfilesDir/optimize_augustus.stderr";
                $stdoutfile = "$otherfilesDir/optimize_augustus.stdout";
	            if($nice){
	                $gb_good_size = `nice grep -c LOCUS $otherfilesDir/genbank.good.gb`;
	            }else{
                    $gb_good_size = `grep -c LOCUS $otherfilesDir/genbank.good.gb`;
	            }
	            $perlCmdString = "";
	            if($nice){
	                $perlCmdString .= "nice ";
	            }
                if($gb_good_size <= 1000){
	                if($nice){
	    	            $perlCmdString .= "perl $string --nice --species=$species --cpus=$CPU --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH $otherfilesDir/genbank.good.gb 1>$stdoutfile 2>$errorfile";
    	            }else{
	    	            $perlCmdString .= "perl $string --species=$species --cpus=$CPU --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH $otherfilesDir/genbank.good.gb 1>$stdoutfile 2>$errorfile";
	                }
                }else{
	                if($nice){
	    	            $perlCmdString .= "perl $string --nice --species=$species --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH --onlytrain=$otherfilesDir/genbank.good.gb.train --cpus=$CPU $otherfilesDir/genbank.good.gb.test 1>$stdoutfile 2>$errorfile";
	                }else{
	    	            $perlCmdString .= "perl $string  --species=$species --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH --onlytrain=$otherfilesDir/genbank.good.gb.train --cpus=$CPU $otherfilesDir/genbank.good.gb.test 1>$stdoutfile 2>$errorfile";
	                }
                }
                print LOG "\# ".(localtime).": optimize AUGUSTUS parameter\n";
                print LOG "$perlCmdString\n\n";
                system("$perlCmdString")==0 or die("failed to execute: $perlCmdString!\n");
                print STDOUT "parameter optimized.\n";
            }
        }

        # train AUGUSTUS for the second time
        if(!uptodate(["$otherfilesDir/genbank.good.gb.train","$otherfilesDir/genbank.good.gb"],["$otherfilesDir/secondetraining.stdout"])){
            print STDOUT "NEXT STEP: second etraining\n";
            $augpath = "$AUGUSTUS_BIN_PATH/etraining";
            $errorfile = "$errorfilesDir/secondetraining.stderr";
            $stdoutfile = "$otherfilesDir/secondetraining.stdout";
            if($nice){
	            $gb_good_size = `nice grep -c LOCUS $otherfilesDir/genbank.good.gb`;
            }else{
                $gb_good_size = `grep -c LOCUS $otherfilesDir/genbank.good.gb`;
            }
            $cmdString = "";
            if($nice){
	            $cmdString .= "nice ";
            }
            if($gb_good_size <= 1000){
                $cmdString .= "$augpath --species=$species --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH $otherfilesDir/genbank.good.gb 1>$stdoutfile 2>$errorfile";
            }else{
                $cmdString .= "$augpath --species=$species --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH $otherfilesDir/genbank.good.gb.train 1>$stdoutfile 2>$errorfile";
            }
            print LOG "\# ".(localtime).": second etraining\n";
            print LOG "$cmdString\n\n";
            system("$cmdString")==0 or die("failed to execute etraining (second time): $!\n");
            print STDOUT "second etraining complete\n";
        }

        # second test
        if(!uptodate(["$otherfilesDir/genbank.good.gb.test","$otherfilesDir/genbank.good.gb"],["$otherfilesDir/secondtest.out"]) || $overwrite){
            $augpath = "$AUGUSTUS_BIN_PATH/augustus";
            print STDOUT "NEXT STEP: second test\n";
            $errorfile = "$errorfilesDir/secondtest.stderr";
            $stdoutfile = "$otherfilesDir/secondtest.stdout";
            if($nice){
                $gb_good_size = `nice grep -c LOCUS $otherfilesDir/genbank.good.gb`;
            }else{
	            $gb_good_size = `grep -c LOCUS $otherfilesDir/genbank.good.gb`;
            }
            $cmdString = "";
            if($nice){
	            $cmdString .= "nice ";
            }
            if($gb_good_size <= 1000){
                $cmdString .= "$augpath --species=$species --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH $otherfilesDir/genbank.good.gb >$stdoutfile 2>$errorfile";
            }else{
                $cmdString .= "$augpath --species=$species --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH $otherfilesDir/genbank.good.gb.test >$stdoutfile 2>$errorfile";
            }
            print LOG "\# ".(localtime).": second test\n";
            print LOG "$cmdString\n\n";
            system("$cmdString")==0 or die("failed to execute augustus (second test): $!\n");
            $target_2 = accuracy_calculator($stdoutfile);
	        print LOG "\# The accuracy after training (after optimize_augustus.pl, no CRF) is $target_2\n";    
            print STDOUT "second test finished.\n";              
        }

        # optional CRF training
        if($crf){
	        if(!uptodate(["$otherfilesDir/genbank.good.gb.test","$otherfilesDir/genbank.good.gb"],["$otherfilesDir/crftraining.stdout"]) || $overwrite){
	            $augpath = "$AUGUSTUS_BIN_PATH/etraining";
            }
	        print STDOUT "NEXT STEP: CRF training\n";
	        $errorfile = "$errorfilesDir/crftraining.stderr";
	        $stdoutfile = "$otherfilesDir/crftraining.stdout";
	        # THIS IS WHERE SIMONE STOPPED
	    
	        if($nice){
	            $gb_good_size = `nice grep -c LOCUS $otherfilesDir/genbank.good.gb`;
            }else{
                $gb_good_size = `grep -c LOCUS $otherfilesDir/genbank.good.gb`;
            }
            $cmdString = "";
            if($nice){
	            $cmdString .= "nice ";
            }
            if($gb_good_size <= 1000){
                $cmdString .= "$augpath --species=$species --CRF=1 --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH $otherfilesDir/genbank.good.gb 1>$stdoutfile 2>$errorfile";
            }else{
                $cmdString .= "$augpath --species=$species --CRF=1 --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH $otherfilesDir/genbank.good.gb.train 1>$stdoutfile 2>$errorfile";
            }
            print LOG "\# ".(localtime).": third etraining: with CRF\n";
            print LOG "$cmdString\n\n";
            system("$cmdString")==0 or die("failed to execute etraining (third time, with CRF): $!\n");
            print STDOUT "second etraining (CRF) complete\n";
        
            # third test
            if(!uptodate(["$otherfilesDir/genbank.good.gb.test","$otherfilesDir/genbank.good.gb"],["$otherfilesDir/thirdtest.out"]) || $overwrite){
                $augpath = "$AUGUSTUS_BIN_PATH/augustus";
                print STDOUT "NEXT STEP: third test\n";
                $errorfile = "$errorfilesDir/thirdtest.stderr";
                $stdoutfile = "$otherfilesDir/thirdtest.stdout";
                if($nice){
                    $gb_good_size = `nice grep -c LOCUS $otherfilesDir/genbank.good.gb`;
                }else{
	                $gb_good_size = `grep -c LOCUS $otherfilesDir/genbank.good.gb`;
                }
                $cmdString = "";
                if($nice){
	                $cmdString .= "nice ";
                }
                if($gb_good_size <= 1000){
                    $cmdString .= "$augpath --species=$species --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH $otherfilesDir/genbank.good.gb >$stdoutfile 2>$errorfile";
                }else{
                    $cmdString .= "$augpath --species=$species --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH $otherfilesDir/genbank.good.gb.test >$stdoutfile 2>$errorfile";
                }
                print LOG "\# ".(localtime).": third test\n";
                print LOG "$cmdString\n\n";
                system("$cmdString")==0 or die("failed to execute augustus (third test): $!\n");
                $target_3 = accuracy_calculator($stdoutfile);
	            print LOG "\# The accuracy after CRF training is $target_2\n";    
                print STDOUT "third test finished.\n";              
            }
            
            # decide on whether to keep CRF parameters
            if($target_2>$target_3){
	            print "CRF performance is worse than HMM performance\n";
	            print LOG "\n\n####### CRF performance is worse than HMM performance #######\n\n\n";
	        }
	        # cp config files
	        print "Copy parameter files to $species*.CRF\n";
	        for(("$species"."_exon_probs.pbl","$species"."_igenic_probs.pbl", "$species"."_intron_probs.pbl")){
	            $cmdString = "cp $AUGUSTUS_CONFIG_PATH/species/$species/$_ cp $AUGUSTUS_CONFIG_PATH/species/$species/$_".".CRF";
	            system($cmdString)==0 or die("failed to execute: copying of CRF parameters $!\n");
	        }
	        # if the accuracy doesn't improve with CRF, overwrite the config files with the HMM parameters from last etraining
	        if($target_2>$target_3){
	            for(("$species"."_exon_probs.pbl","$species"."_igenic_probs.pbl", "$species"."_intron_probs.pbl")){
	                $cmdString = "rm $AUGUSTUS_CONFIG_PATH/species/$species/$_";
                    system("$cmdString")==0 or die("failed to execute: removing CRF parameter file $_!\n");
                    $cmdString = "cp $AUGUSTUS_CONFIG_PATH/species/$species/$_".".HMM $AUGUSTUS_CONFIG_PATH/species/$species/$_";
                    system("$cmdString")==0 or die("faled to execute: copying withoutCRF parameter file to parameter file in use!\n");
		        }
	        }	          
        }
    }
  
    # copy species files to working directory
    if(! -d "$parameterDir/$species"){
        print STDOUT "NEXT STEP: copy optimized parameters to working directory\n";
        $cmdString = "";
        if($nice){
	        $cmdString .= "nice ";
        }
        $cmdString .= "cp -r $AUGUSTUS_CONFIG_PATH/species/$species $parameterDir";
        print LOG "\# ".(localtime).": copy optimized parameters to working directory\n";
        print LOG "$cmdString\n\n";
        system("$cmdString")==0 or die("failed to execute: $cmdString!\n");
        print STDOUT "parameter files copied.\n";
    }
}



         ####################### run AUGUSTUS  #########################
# run AUGUSTUS for given species with given options
sub augustus{
  my $scriptpath;
  $augpath = "$AUGUSTUS_BIN_PATH/augustus";
  if(defined($augustus_scripts_path)){
    $scriptpath = $AUGUSTUS_SCRIPTS_PATH;
  }elsif(defined($AUGUSTUS_BIN_PATH)){
    $scriptpath = "$AUGUSTUS_BIN_PATH/../scripts";
  }else{
    $scriptpath = "$AUGUSTUS_CONFIG_PATH/../scripts";
  }
  my $extrinsic;
  my @genome_files;
  my $pm;
  if(defined($extrinsicCfgFile) && -f $extrinsicCfgFile){
    $extrinsic = $extrinsicCfgFile;
    print STDOUT "Use assigned $extrinsicCfgFile as extrinsic file\n";
  }elsif($useexisting || ! -f "$AUGUSTUS_CONFIG_PATH/species/$species/extrinsic.$species.cfg"){
    $extrinsic = "$AUGUSTUS_CONFIG_PATH/extrinsic/extrinsic.M.RM.E.W.cfg";
  }else{
    $extrinsic = "$AUGUSTUS_CONFIG_PATH/species/$species/extrinsic.$species.cfg";
  }
 
  if(!uptodate([$extrinsic,$hintsfile, $genome],["$otherfilesDir/augustus.gff"])  || $overwrite){
    # split genome file in smaller parts and use multiple parallel runs of augustus for prediction
    if($CPU > 1){
      print STDOUT "NEXT STEP: split genome file in smaller parts\n";
      $string = find("splitMfasta.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);    
      $errorfile = "$errorfilesDir/splitMfasta.stderr";
      my $minsize = floor($genome_length / $CPU);
      $perlCmdString = "";
      if($nice){
	  $perlCmdString .= "nice ";
      }
      $perlCmdString .= "perl $string $genome --outputpath=$otherfilesDir --minsize=$minsize 2>$errorfile";
      system("$perlCmdString")==0 or die("failed to execute: $perlCmdString!\n");
      if($nice){
          @genome_files = `nice find $otherfilesDir -name "genome.split.*"`;
      }else{
	  @genome_files = `find $otherfilesDir -name "genome.split.*"`;
      }
      print LOG "\# ".(localtime).": split genome file in ".scalar(@genome_files)." parts\n";
      print LOG "$perlCmdString\n\n";
      print STDOUT "split genome file in ".scalar(@genome_files)." parts complete.\n";
    }else{
      push(@genome_files, $genome);
    }
    @genome_files = sort {lc($a) cmp lc($b)} @genome_files;
    $pm = new Parallel::ForkManager($CPU);
    print STDOUT "NEXT STEP: run AUGUSTUS\n";
    for(my $i = 0; $i < scalar(@genome_files); $i++){
      chomp($genome_files[$i]);
      my $pid = $pm->start and next;
      my $idx = $i + 1;
      $errorfile = "$errorfilesDir/augustus.$idx.stderr";
      $stdoutfile = "$otherfilesDir/augustus.$idx.gff";
      $cmdString = "";
      if($nice){
	  $cmdString .= "nice ";
      }
      $cmdString .= "$augpath --species=$species --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH --extrinsicCfgFile=$extrinsic --alternatives-from-evidence=$alternatives_from_evidence --hintsfile=$hintsfile --UTR=$UTR";
      if(defined($optCfgFile)){
        $cmdString .= " --optCfgFile=$optCfgFile"; 
      }
      if($soft_mask){
        $cmdString .= " --softmasking=1";
      }
      if(defined($augustus_args)){
	  $cmdString .= " $augustus_args";
      }
      $cmdString .= " $genome_files[$i] 1>$stdoutfile 2>$errorfile";
      print LOG "\# ".(localtime).": run AUGUSTUS for file $genome_files[$idx-1]\n";
      print LOG "$cmdString\n\n";
      system("$cmdString")==0 or die("failed to execute augustus (prediction step): $!\n");
      $pm->finish;

    }
    $pm->wait_all_children;
    print STDOUT "AUGUSTUS complete.\n";
    # join prediction files to one file
    if($CPU > 1){
      print STDOUT "NEXT STEP: concatenate and join AUGUSTUS output files\n";
      $string = find("join_aug_pred.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
      my $cat_file = "$otherfilesDir/augustus.tmp.gff";
      for(my $idx = 1; $idx <= scalar(@genome_files); $idx++){
	  $cmdString = "";
	  if($nice){
	      $cmdString .= "nice ";
	  }
	  $cmdString .= "cat $otherfilesDir/augustus.$idx.gff >> $cat_file";
	  system("$cmdString")==0 or die("failed to execute cat augustus.idx.gff files: $!\n");
      }
      $cmdString = "";
      if($nice){
	  $cmdString .= "nice ";
      }
      $cmdString .= "perl $string <$cat_file >$otherfilesDir/augustus.gff";
      print LOG "\# ".(localtime).": concatenate and join AUGUSTUS output files\n";
      print LOG "$cmdString\n\n";
      system("$cmdString")==0 or die("failed to execute join_aug_pred.pl: $!\n");
      print STDOUT "join AUGUSTUS predictions complete.\n";
      for(my $idx = 1; $idx <= scalar(@genome_files); $idx++){
	  unlink("$otherfilesDir/augustus.$idx.gff");
	  print LOG "Deleting $otherfilesDir/augustus.$idx.gff\n";
      }
      unlink("$otherfilesDir/augustus.tmp.gff");
      print LOG "Deleting $otherfilesDir/augustus.tmp.gff\n";
    }else{
	$cmdString = "";
	if($nice){
	    $cmdString .= "nice ";
	}
	$cmdString .= "mv $otherfilesDir/augustus.1.gff $otherfilesDir/augustus.gff";
	print LOG "\# ".(localtime).": rename AUGUSTUS file\n";
	print LOG "$cmdString\n\n";
	system("$cmdString")==0 or die("failed to execute mv for renaming AUGUSTUS file: $!\n");
    }
  }
}

# call getAnnoFasta.pl
sub getAnnoFasta{
  my $AUG_pred = shift;
  @_ = split(/\//, $AUG_pred);
  my $name_base = substr($_[-1],0,-4);
  my $string = find("getAnnoFasta.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
  my $errorfile = "$errorfilesDir/getAnnoFasta.$name_base.stderr";
  $perlCmdString = "";
  if($nice){
      $perlCmdString .= "nice ";
  }
  my $perlCmdString .= "perl $string $AUG_pred 2>$errorfile";
  print LOG "\# ".(localtime).": make a fasta file with protein sequences for $AUG_pred\n";
  print LOG "$perlCmdString\n\n";
  system("$perlCmdString")==0 or die("failed to execute getAnnoFasta.pl: $!\n");
}

# make gtf file
sub make_gtf{
  my $AUG_pred = shift;
  @_ = split(/\//, $AUG_pred);
  my $name_base = substr($_[-1],0,-4);
  my $gtf_file = substr($AUG_pred,0,-4).".gtf";
  my $errorfile = "$errorfilesDir/gtf2gff.$name_base.gtf.stderr";
  my $perlstring = find("gtf2gff.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
  $cmdString = "";
  if($nice){
      $cmdString .= "nice ";
  }
  my $cmdString .= "cat $AUG_pred | perl -ne 'if(m/\\tAUGUSTUS\\t/){print \$_;}' | perl $perlstring --printExon --out=$gtf_file 2>$errorfile";
  print "$cmdString\n\n";
  print LOG "\# ".(localtime).": make a gtf file from $AUG_pred\n";
  print LOG "$cmdString\n\n";
  system("$cmdString")==0 or die("failed to execute  gtf2gff.pl: $!\n");
  if($gff3){
    my $gff3_file = substr($AUG_pred,0,-4).".gff3";
    my $errorfile = "$errorfilesDir/gtf2gff.$name_base.gff3.stderr";
    my $perlstring = find("gtf2gff.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
    $cmdString = "";
    if($nice){
	$cmdString .= "nice ";
    }
    my $cmdString .= "cat $AUG_pred | perl -ne 'if(m/\\tAUGUSTUS\\t/){print \$_;}' | perl $perlstring --printExon -gff3 --out=$gff3_file 2>$errorfile";
    print "$cmdString\n\n";
    print LOG "\# ".(localtime).": make a gff3 file from $AUG_pred\n";
    print LOG "$cmdString\n\n";
    system("$cmdString")==0 or die("failed to execute gtf2gff.pl: $!\n");
  }
}   
  

         ####################### delete all zero sized files #########################
# delete empty files
sub clean_up{
    print STDOUT "NEXT STEP: delete empty files\n";
    if($nice){
	@files = `nice find $otherfilesDir -empty`;
    }else{
        @files = `find $otherfilesDir -empty`;
    }
  print LOG "\# ".(localtime).": delete empty files\n";
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
# checks for GeneMark-ET: perl modules: YAML, Hash::Merge, Logger::Simple, Parallel::ForkManager
# checks for braker: perl modules: Scalar::Util::Numeric
sub check_upfront{ # see autoAug.pl

  # check whether required perl modules are installed
  my $pmodule;
  my @module_list = (
    "YAML",
    "Hash::Merge",
    "Logger::Simple",
    "Parallel::ForkManager",
    "Scalar::Util::Numeric"
  );

  foreach my $module (@module_list){  
    $pmodule = check_install(module => $module);
    if(!$pmodule){
      print STDOUT "WARNING: Perl module '$module' is required but not installed yet.\n";
    }
  }

  # check whether environmental variables are set or other path variables are defined
  if(!$ENV{'AUGUSTUS_CONFIG_PATH'} && !defined($augustus_cfg_path)){ # see autoAug.pl
    print STDERR "ERROR: The environment variable AUGUSTUS_CONFIG_PATH is not defined. Please export an environment variable for AUGUSTUS or use --AUGUSTUS_CONFIG_PATH=path/to/augustus/config .\n"; # see autoAug.pl
    exit(1);
  }

  if(!defined($augustus_bin_path)){
    print STDOUT "WARNING: The command line option --AUGUSTUS_BIN_PATH was not used. This is ok if binaries of augustus reside in ../bin relative to AUGUSTUS_CONFIG_PATH. Otherwise, please specify --AUGUSTUS_BIN_PATH\n";
  }
  
  if(!defined($augustus_scripts_path)){
    print STDOUT "WARNING: The command line option --AUGUSTUS_SCRIPTS_PATH was not used. This is ok if scripts of augustus reside in ../scripts relative to AUGUSTUS_CONFIG_PATH or AUGUSTUS_BIN_PATH (will first look for them in AUGUSTUS_BIN_PATH). Otherwise, please specify --AUGUSTUS_SCRIPTS_PATH\n";
  }
  
  if(!$ENV{'GENEMARK_PATH'} && !defined($GMET_path)){
    print STDERR "ERROR: The environment variable GENEMARK_PATH to the 'gmes_petap.pl' script is not defined. Please export an environment variable or use --GENEMARK_PATH=path/to/gmes_petap.pl.\n"; 
    exit(1);
  }

  if(!$ENV{'ALIGNMENT_TOOL_PATH'} && !defined($ALIGNMENT_TOOL_PATH) && @prot_seq_files){
      print STDERR "ERROR: The environment variable  ALIGNMENT_TOOL_PATH to either the exectuable of GenomeThreader, Exonerate or Spaln is not defined. Please expoert an environment variable or use --ALIGNMENT_TOOL_PATH=path/to/aligner.\n";
  }

  if(!$ENV{'BAMTOOLS_PATH'} && !defined($bamtools_path)){ # see autoAug.pl
    print STDERR "ERROR: The environment variable BAMTOOLS_PATH is not defined. Please export an environment variable for bamtools or use --BAMTOOLS_PATH=path/to/bamtools.\n"; # see autoAug.pl
    exit(1);
  }

  # check for augustus executable
  $augpath = "$AUGUSTUS_BIN_PATH/augustus";
  if(system("$augpath > /dev/null 2> /dev/null") != 0){                   # see autoAug.pl
    if(! -f $augpath){                                                    # see autoAug.pl
      print STDERR "ERROR: augustus executable not found at $augpath.\n"; # see autoAug.pl
    }else{
      print STDERR "ERROR: $augpath not executable on this machine.\n";   # see autoAug.pl
    }
    exit(1);
  }
  
  # check whether bamtools is installed
  if(system("which $BAMTOOLS_BIN_PATH/bamtools > /dev/null") != 0){
    print STDERR "Error: bamtools not installed. Please install it first.\n";
    exit (1);
  }

  # check for etraining executable
  my $etrainpath;
  $etrainpath = "$AUGUSTUS_BIN_PATH/etraining";
  if(system("$etrainpath > /dev/null 2> /dev/null") != 0){                   
    if(! -f $etrainpath){                                                    
      print STDERR "ERROR: etraining executable not found at $etrainpath.\n";
    }else{
      print STDERR "ERROR: $etrainpath not executable on this machine.\n";
    }
    exit(1);
  }

  # check for alignment executable and in case of SPALN for environment variables
  my $prot_aligner;
  if(@prot_seq_files){
      if($prg eq 'gth'){
	  $prot_aligner = "$ALIGNMENT_TOOL_PATH/gth";
	  if(! -f $prot_aligner){
	      print STDERR "ERROR: GenomeThreader executable not found at $prot_aligner.\n";
	      exit(1);
	   }elsif(! -x $prot_aligner){
	      print STDERR "ERROR: $prot_aligner not executable on this machine.\n";
	      exit(1);
	  }
      }elsif($prg eq 'spaln'){
	  $prot_aligner = "$ALIGNMENT_TOOL_PATH/spaln";
	  if(! -f $prot_aligner){
	      print STDERR "ERROR: Spaln executable not found at $prot_aligner.\n";
	      exit(1);
	  }elsif(! -x $prot_aligner){
	      print STDERR "ERROR: $prot_aligner not executable on this machine.\n";
	      exit(1);
	  }
	  # check whether spaln environment variables are configured
	  if(!$ENV{'ALN_DBS'} or !$ENV{'ALN_TAB'}){
	      if(!$ENV{'ALN_DBS'}){
		  print STDERR "ERROR: The environment variable ALN_DBS for spaln is not defined. Please export an\
 environment variable with:' export ALN_DBS=/path/to/spaln/seqdb'\n";
	      }
	      if(!$ENV{'ALN_TAB'}){
		  print STDERR "ERROR: The environment variable ALN_TAB for spaln is not defined. Please export an\
 environment variable with:' export ALN_TAB=/path/to/spaln/table'\n";
	      }
	      exit(1);
	  }
      }elsif($prg eq 'exonerate'){
	  $prot_aligner = "$ALIGNMENT_TOOL_PATH/exonerate";
	  if(! -f $prot_aligner){
	      print STDERR "ERROR: Exonerate executable not found at $prot_aligner.\n";
	      exit(1);
	  }elsif(! -x $prot_aligner){
	      print STDERR "ERROR: $prot_aligner not executable on this machine.\n";
	      exit(1);
	  }
	  
      }
  }
  
  
  # check whether the necessary perl scripts exist and can be found
  find("gff2gbSmallDNA.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
  find("filterGenemark.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
  find("filterIntronsFindStrand.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
  find("new_species.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
  find("filterGenesIn_mRNAname.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
  find("join_mult_hints.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
  find("randomSplit.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
  find("optimize_augustus.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
  find("splitMfasta.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
  find("join_aug_pred.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
  find("getAnnoFasta.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
  find("gtf2gff.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
  find("startAlign.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
  find("align2hints.pl", $AUGUSTUS_BIN_PATH, $AUGUSTUS_SCRIPTS_PATH, $AUGUSTUS_CONFIG_PATH);
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
    my $isAllowed = 0;
    foreach(@allowedHints){
	if($gff_line[2] eq $_){
	    $isAllowed = 1;
	}
    }
    if($isAllowed != 1){
	print STDERR "ERROR: File $gfffile contains hints of a feature type that is currently not supported by BRAKER\n";
	print STDERR "Currently allowed hint types:\n";
	foreach(@allowedHints){
	    print STDERR $_."\n";
	}
	close(GFF) or die("Could not close gff file $gfffile!\n");
	exit(1);
    }
  }
  close(GFF) or die("Could not close gff file $gfffile!\n");
  print STDOUT "gff format check complete.\n";
}

# check whether all options are set correctly
sub check_options{
  print STDOUT "NEXT STEP: check options\n";
  if($alternatives_from_evidence ne "true" && $alternatives_from_evidence ne "false"){
    print STDERR "ERROR: \"$alternatives_from_evidence\" is not a valid option for --alternatives-from-evidence. Please use either 'true' or 'false'.\n";
    exit(1);
  } 

  if($UTR ne "on" && $UTR ne "off"){
    print STDERR "ERROR: \"$UTR\" is not a valid option for --UTR. Please use either 'on' or 'off'.\n";
    exit(1);
  }

  my $operatingSystem = "$^O";
  my $cpus_available = 1;
  if($operatingSystem eq "linux"){
      $cpus_available = `nproc`;
  }else{ # Mac OS X
      $cpus_available = `sysctl -n hw.ncpu`;
  }

  if($cpus_available < $CPU){
    print STDOUT "WARNING: Your system does not have $CPU cores available, only $cpus_available. Braker will use the $cpus_available available instead of the chosen $CPU.\n";
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
  my $stdStr = "This may later on cause problems! The pipeline will create a new file without spaces or \"|\" characters and a header.map file to look up the old and new headers. This message will be suppressed from now on!\n";
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
        # replace | and whitespaces by _
	my $oldHeader = $_;
	$_ =~ s/\s/_/g;
	$_ =~ s/\|/_/g;
	print OUTPUT "$_\n";
	print MAP "$_\t$oldHeader\n";
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

# check bam headers
sub check_bam_headers{
  
  my $bamFile = shift;
  my $someThingWrongWithHeader = 0; # see autoAug.pl
  my $spaces = 0;                   # see autoAug.pl
  my $orSign = 0;                   # see autoAug.pl
  my %map_hash;
  my $mapFile = "$otherfilesDir/bam_header.map";
  my $stdStr = "This may later on cause problems! The pipeline will create a new file without spaces or \"|\" characters and a bam_header.map file to look up the old and new headers, if samtools is working on your system. This message will be suppressed from now on!\n";

  @_ = split(/\//, $bamFile); 
  @_ = split(/\./, $_[-1]);
  my $samHeaderFile = "$otherfilesDir/".$_[0]."_header.sam";
  my $samHeaderFile_new = "$otherfilesDir/".$_[0]."_new_header.sam";
  
  if(!uptodate([$bamFile],["$otherfilesDir/$bamFile"]) || $overwrite){
    # extract header information 
      print STDOUT "NEXT STEP: create SAM header file $samHeaderFile.\n";
      $cmdString = "";
      if($nice){
	  $cmdString .= "nice ";
      }
      $cmdString .= "$BAMTOOLS_BIN_PATH/bamtools header -in $bamFile > $samHeaderFile";
      print LOG "\# ".(localtime).": create header file $samHeaderFile\n";
      print LOG "$cmdString\n\n";
      system("$cmdString")==0 or die("failed to execute: $cmdString!\n");
      print STDOUT "SAM file $samHeaderFile complete.\n";

      print STDOUT "NEXT STEP: check BAM headers\n";
      open(SAM, "<", $samHeaderFile) or die("Could not open SAM file $samHeaderFile!\n");
      open(OUTPUT, ">", "$samHeaderFile_new") or die("Could not open SAM file $samHeaderFile_new!\n");
      open(MAP, ">", $mapFile) or die("Could not open map file $mapFile.\n");
      while(<SAM>){
	  chomp;
	  # only check sequence entries
	  if($_ =~ m/^\@SQ/){
	      my @seq_line = split(/\t/, $_);
	      my $seq_end = $seq_line[-1];

	      @seq_line = split(/\:/, $seq_line[1]);
	      my $old_name = $seq_line[1];
	      my $new_name = $old_name;
	      # remove whitespaces, if necessary
	      @seq_line = split(/\s/, $seq_line[1]);

	      if(scalar(@seq_line) > 1){
		  if($spaces == 0){
		      print STDOUT "WARNING: Detected whitespace in BAM header of file $bamFile. ".$stdStr; # see autoAug.pl
		      print STDOUT "Replacing whitespaces by underscores in Bam headers.\n";
		      $spaces++;        # see autoAug.pl
		  }
	      }
	      $new_name =~ s/\s/_/g; # removing whitespaces (if any)
	      @seq_line = split(/\|/, $old_name);
	      if(scalar(@seq_line) > 1){
		  if($orSign == 0){   # see autoAug.pl
		      print STDOUT "WARNING: Detected | in header of file $bamFile. ".$stdStr; # see autoAug.pl
		      print STDOUT "Replacing | by underscores in Bam headers.\n";
		      $orSign++;        # see autoAug.pl
		  }
	      }
	      $new_name =~ s/\|/_/g; # replace or signs by underscores (if any)
	      $map_hash{$old_name} = $new_name;
	      $seq_line[0] = "\@SQ\tSN:$new_name\t$seq_end";
	      if($seq_line[0] !~ m/[>a-zA-Z0-9]/){
		  if($someThingWrongWithHeader==0){   # see autoAug.pl
		      print STDOUT "WARNING: BAM headers in file $bamFile seem to contain non-letter and non-number characters. That means they may contain some kind of special character. ".$stdStr; # see autoAug.pl
		      $someThingWrongWithHeader++;      # see autoAug.pl
		  }
	      }

	      print OUTPUT "$seq_line[0]\n";   # see simplifyFastaHeaders.pl
	      print MAP "$map_hash{$old_name}\t$old_name\n";  # see simplifyFastaHeaders.pl
	  }elsif(eof){
	      print OUTPUT "$_";
	  }else{
	      print OUTPUT "$_\n";
	  }
    }
    close(SAM) or die("Could not close header file $samHeaderFile!\n");
    close(OUTPUT) or die("Could not close output SAM file $samHeaderFile_new!\n");
    close(MAP) or die("Could not close map file $mapFile!\n");
    print STDOUT "headers check for BAM file $bamFile complete.\n";
    print STDOUT "Deleting SAM header file $samHeaderFile (will not be needed from here on)\n";
    unlink($samHeaderFile);
    # something wrong with header part
    if($spaces != 0 || $orSign != 0 || $someThingWrongWithHeader != 0){
      # no samtools installed. stop here
      if(system("which samtools > /dev/null") != 0){
        print STDOUT "'samtools' not installed. BAM file cannot be fixed automatically.\n";
        print STDERR "BAM file $bamFile contains spaces, \"|\" or some other kind of special characters.\n";
        exit(1);
      # samtools installed. try to correct BAM file
      }else{
        if(!$ENV{'SAMTOOLS_PATH'} && !defined($SAMTOOLS_PATH_OP)){
          print STDOUT "WARNING: The environment variable SAMTOOLS_PATH is not defined. Please export an environment variable for samtools or use --SAMTOOLS_PATH=path/to/samtools.\n"; # see autoAug.pl
          print STDOUT "The program will try to use 'samtools' to start samtools, which may not work on your system.\n"; # see autoAug.pl
          $SAMTOOLS_PATH = "samtools";
        }
        print STDOUT "NEXT STEP: convert BAM file to SAM format.\n";
        my $samFile = "$otherfilesDir/".$_[0].".sam";
        my $samFile_new = "$otherfilesDir/".$_[0]."_new.sam";
	$cmdString = "";
	if($nice){
	    $cmdString .= "nice ";
	}
        $cmdString .= "$SAMTOOLS_PATH view $bamFile > $samFile";
        print LOG "\# ".(localtime).": convert BAM to SAM file $samFile\n";
        print LOG "$cmdString\n\n";
        system("$cmdString")==0 or die("failed to execute: $cmdString!\n");
        print STDOUT "SAM file $samFile created.\n";
        open(SAM, "<", $samFile) or die("Could not open SAM file $samFile!\n");
        open(OUTPUT, ">", "$samFile_new") or die("Could not open SAM file $samFile_new!\n");
        while(<SAM>){
          chomp;
          my @line = split(/\t/, $_);
          $line[2] = $map_hash{$line[2]};
          if(eof){
            print OUTPUT join("\t", @line);
          }else{
            print OUTPUT join("\t", @line)."\n";
          }
        }
        close(SAM) or die("Could not close SAM file $samFile!\n");
        close(OUTPUT) or die("Could not close output SAM file $samFile_new!\n");
        print STDOUT "NEXT STEP: concatenate new header and SAM file.\n";
	$cmdString = "";
	if($nice){
	    $cmdString .= "nice ";
	}
        $cmdString .= "cat $samHeaderFile_new $samFile_new > $samFile";
        print LOG "\# ".(localtime).": concatenate new header and SAM file\n";
        print LOG "$cmdString\n\n";
        system("$cmdString")==0 or die("failed to execute: $cmdString!\n");
        print STDOUT "new SAM file created.\n";
        unlink($samFile_new);

        print STDOUT "NEXT STEP: convert new SAM file to BAM format.\n";
	$cmdString = "";
	if($nice){
	    $cmdString .= "nice ";
	}
        $cmdString = "$SAMTOOLS_PATH view -bSh $samFile > $otherfilesDir/".$_[0].".bam";
        print LOG "\# ".(localtime).": convert new SAM file to BAM format\n";
        print LOG "$cmdString\n\n";
        system("$cmdString")==0 or die("failed to execute: $cmdString!\n");
        print STDOUT "new BAM file created.\n";
        unlink($samFile);
        $bamFile = "$otherfilesDir/".$_[0].".bam";
      }
    }
    unlink($samHeaderFile_new);
  } 
  return $bamFile;
} 

# calculate the result of testing AUGUSTUS on genbank files in a single number 
sub accuracy_calculator{
    my $aug_out=shift;
    my ($nu_sen, $nu_sp, $ex_sen, $ex_sp, $gen_sen, $gen_sp);
    open(AUGOUT, "$aug_out") or die ("Could not open $aug_out!\n");
    while(<AUGOUT>){
        if(/^nucleotide level\s*\|\s*(\S+)\s*\|\s*(\S+)/){
            $nu_sen=$1;
            $nu_sp=$2;
        }
        if(/^exon level\s*\|.*\|.*\|.*\|.*\|.*\|\s*(\S+)\s*\|\s*(\S+)/){
            $ex_sen=$1;
            $ex_sp=$2;
        }
        if(/^gene level\s*\|.*\|.*\|.*\|.*\|.*\|\s*(\S+)\s*\|\s*(\S+)/){
            $gen_sen=$1;
            $gen_sp=$2;
        }
    }
    my $target=(3*$nu_sen+2*$nu_sp+4*$ex_sen+3*$ex_sp+2*$gen_sen+1*$gen_sp)/15;
    return $target;
}



