package helpMod;

        ###########################
        # Some help sub functions #
        ###########################

use Exporter 'import';
@EXPORT_OK = qw(find tildeConvert checkFile formatDetector relToAbs setParInConfig uptodate);

use strict;
use Cwd;
use File::Spec::Functions qw(rel2abs);
use File::Basename qw(dirname);

###################################################################################################
# search a script under $AUGUSTUS_CONFIG_PATH and the directory where this script placed  in turn #
# return the name with absolute path of this script                                               #
# usage: $absNameOfScript=find("script")                                                          #
###################################################################################################

sub find{
  my $script=shift;           # script to find
  my $exist;                  # boolean variable, to remark if $script is found
  my $string;                 # the absolute name for $script

  my $AUGUSTUS_CONFIG_PATH = $ENV{'AUGUSTUS_CONFIG_PATH'};      # the environment varialbe AUGUSTUS_CONFIG_PATH
  my $path_1="$AUGUSTUS_CONFIG_PATH/../scripts";                # first searching path of $script
  my $path_2=dirname(rel2abs($0));                              # second searching path of $script

  foreach(("$path_1/$script", "$path_2/$script")){
    if(-f $_){
      $exist=1;
      $string=$_;
      last;
    }
  }
  if($exist){
    #print "Found script $string.\n";
    return $string;                    # return name with absolute path 
  }
  else{
    # if not found, output error
    die("Error: found neither $path_1/$script nor $path_2/$script!\nPlease Check the environment variable AUGUSTUS_CONFIG_PATH or install AUGUSTUS again!\n");
  }
}

######################################################################
# convert a file name which begins with ~ to name with absolute path #
######################################################################

sub tildeConvert{
    my $file=shift;
    if($file =~ /^~/){
	my $HOME = $ENV{'HOME'};     
	$file=~s/~/$HOME/;  # replace ~ with home directory
	
    }
    return $file;
}

##################################################################
# check if $file exists and replace $file with its absolute path #   
##################################################################

sub checkFile{
    my $file=shift;      # file which to be checked
    my $type=shift;      # type of file, used by error outputting if the file doesn't exist
    my $usage=shift;     # usage to be outputted if the file doesn't exist

    die("Error: missing $type file!\n$usage") if(!$file);
    
    # overwrite $file with absolute path
    $file=tildeConvert($file);
    $file=rel2abs($file);  # overwrite $file with absolute path
    if(!(-f $file)){
        die("Error: $type file $file not found!\n");
    }
    return $file;
}

##################################################
# detect if a file has gff or gb or fasta format #
##################################################

sub formatDetector{
    my $file=shift;   # file to be detected
    my $testLines = 1000; # read at most this many lines for testing
    my $i;
    my @helpArray_gff;
    #
    # check if file has GENBANK format
    #
    open(DFILE, $file) or die ("Could not open $file!\n");
    $i=0;
    my $haveLOCUS=0;
    my $haveSource=0;
    my $haveOrigin=0;
    my $haveTermSymb=0;
    while(defined(my $line=<DFILE>) && $i<$testLines){
	$i++;
	$haveLOCUS++ if($i==1 && $line=~ /^LOCUS/);
	$haveSource++ if($line=~ / +source +/i);
	$haveOrigin++ if($line=~ /^ORIGIN/);
	$haveTermSymb++ if($line=~ /\/\//);
    }
    close(DFILE);
    if ((($haveLOCUS>0) + ($haveSource>0) + ($haveOrigin>0) + ($haveTermSymb>0)) > 1){
	print STDERR "$file appears to be in corrupt Genbank format. 'LOCUS' missing\n" if (!$haveLOCUS);
	print STDERR "$file appears to be in corrupt Genbank format. ' source ' line missing\n" if (!$haveSource);
	print STDERR "$file appears to be in corrupt Genbank format. 'ORIGIN' missing\n" if (!$haveOrigin);
	print STDERR "$file appears to be in corrupt Genbank format. '//' missing\n" if (!$haveTermSymb);
	return "gb";
    }
    #
    # check if file has GFF format
    #
    open(DFILE, $file) or die ("Could not open $file!\n");
    $i=0;
    my $badGFFlines=0;
    my $goodGFFlines=0;
    while(defined(my $line=<DFILE>) && $i<$testLines){
	$i++;
        # if not genbank format and the row not a possible comment in gff format
        if(!($line=~/^#/) && !($line=~/^\s*$/)){
             @helpArray_gff=split(/\t/, $line);
	     if($#helpArray_gff<7){ 
               # each non-comment row should contain at least 7 tabulators (end with new line???)
		 $badGFFlines++;
    	     } else {
		 $goodGFFlines++;
	     }
	}
    }
    close(DFILE);
    if ($goodGFFlines > 0){
	if($badGFFlines > 0){
	    print STDERR "$file appears to be in corrupt GFF format.\n";
	    return "";
	} else {
	    return "gff";
	}
    }
    #
    # check if file has FASTA format and whether it is DNA or protein
    #
    open(DFILE, $file) or die ("Could not open $file!\n");
    $i=0;
    my $greaterLines=0;
    my $concatseq = "";
    while(defined(my $line=<DFILE>) && $i<$testLines){
	$i++;
	chomp $line;
	if ($line =~ /^>/){
	    $greaterLines++;
	} else {
	    $concatseq .= $line;
	}
    }
    close(DFILE);
    if ($greaterLines > 0){
	$concatseq = uc($concatseq);
	my $len = length($concatseq);
	my $protchar = ($concatseq =~ tr/ACDEFGHIKLMNPQRSTVWYX//);
	my $dnachar = ($concatseq =~ tr/ACGTN//);
	if ($protchar > 0.8 * $len){
	    return "fasta-prot";
	} elsif ($dnachar > 0.8 * $len){
	    return "fasta-dna";
	} else {
	    print STDERR "$file appears to be in corrupt FASTA format.\n";
	    return "";
	}
    }
    return "";
}

##########################################
# convert relative path to absolute path #
##########################################

sub relToAbs{
    my $name=shift;
    $name=tildeConvert($name);         # overwrite working directory
    return rel2abs($name);             # with absolute path
}


##########################################
# change a parameter in a config file    #
# assume the format                      #
# parName    value   # comment           #
##########################################

sub setParInConfig{
    my $configFileName = shift;
    my $parName = shift;
    my $value = shift;
    open(CFGFILE, "+<$configFileName") or die ("Could not read config file $configFileName\n");
    my @lines = <CFGFILE>;
    foreach my $line (@lines){
	$line =~ s/(\s*$parName +)(\S+?)(\s|\#|$)/$1$value$3/;
    }
    seek(CFGFILE, 0,0);
    print CFGFILE @lines or die ("Could not write $configFileName");
    truncate(CFGFILE, tell(CFGFILE));
    close(CFGFILE);
}

##################################################
# uptodate
# check whether output files are up to date with respect to input files
# all output files must exist and not be older than any input file
##################################################

sub uptodate {
    my $input = shift;    # reference to list of input file names
    my $output = shift;   # reference to list of output file names
    my $earliestOutMtime; # earliest modification time of an output file
    my $latestInMtime;    # latest modification time an any input file
    my @stat;             # holds info about file

    return 1 if (@{$output} == 0); # no output is always up to date
    # check whether all output files exist
    foreach my $of (@{$output}){
	return 0 if (! -f $of);
	@stat = stat($of);
	$earliestOutMtime = $stat[9] if (!defined($earliestOutMtime) || $stat[9] < $earliestOutMtime);
    }
    return 1 if (@{$input} == 0); # no input is always older than output files
    # check existence and times of input files
    foreach my $if (@{$input}){
	if (! -f $if){ # ignore if input file does not exist
	    print STDERR "Warning: $if missing.\n";# TODO, remove or correct this later
	}
	@stat = stat($if);
	$latestInMtime = $stat[9] if (!defined($latestInMtime) || $stat[9] > $latestInMtime);
    }
    return ($latestInMtime <= $earliestOutMtime);
}

1;
