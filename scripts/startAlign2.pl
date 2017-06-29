#!/usr/bin/env perl

####################################################################################################
#                                                                                                  #
# startExonerate.pl - run exonerate on partial sequences and correspondig protein file             #
#                                                                                                  #
# Author: Simone Lange                                                                             #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# Release date: October 28th 2015                                                                  #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
#                                                                                                  #
#################################################################################################### 


# ----------------------------------------------------------------
# | first outline                   | Simone Lange   |28.10.2015 |
# | changed default from region to  |                |01.12.2015 |
# | whole sequence (noreg to reg)   |                |           |
# | added "wholeProt option" (just  |                |           |
# | leave out list and pos files)   |                |           |
# ----------------------------------------------------------------
 
use Getopt::Long;
use Cwd;
use File::Path qw(make_path);
use File::Spec::Functions qw(rel2abs);
use Parallel::ForkManager;

use strict;
use warnings;

my $usage = <<'ENDUSAGE';

startAlign.pl  split genome file in single sequences or sequence parts and protein file according 
               to contigIDs then run alignment programme exonerate, spaln or genomeThreader (gth) 
               for each contigID or each sequence. When no list and/or pos file(s) are/is assigned
               the programme will use the whole protein file.

SYNOPSIS

startAlign.pl [OPTIONS] --genome=genome.fa --prot=db.fa  --prg=gth|exonerate|spaln

  --genome=genome.fa          fasta file with DNA sequences
  --prot=db.fa                fasta file with protein sequences 
  
    
OPTIONS

    --help                       Print this help message
    --CPU=n                      Specifies the maximum number of CPUs that can be used during
    --dir=path/to/dir            Set path to working directory. In the working directory results
                                 and temporary files are stored.
    --list=BLAST.hit.list        Contains contig and protein ID. Format: contigID proteinID
    --log=startAlign.log         Log file
    --maxintronlen=n             Exonerate option: Alignments with longer gaps are discarded (default 30000).
    --reg                        Use region parts and not whole sequences.
    --offset=n                   This many bp are added before and after cutout coordinates.
    --prg=s                      Alignment programme to call. Valid options are 'gth', 'spaln' or 'exonerate'.
    --pos=dna.pos                Contains information on contigs and genome sequence. Format
                                 contigID nr_of_prots_mapped start end strand chrID


DESCRIPTION
      
  Example:

    startAlign.pl [OPTIONS] --genome=genome.fa --prot=db.fa --list=BLAST.hit.list --pos=dna.pos --prg=gth

ENDUSAGE

my $alignDir;               # directory for alignment output files
my $errorfile;              # error file name
my $cmdString;              # to store shell commands
my %contigIDs;              # hash for contig IDs
my $counterW = 0;           # print message 'add amino ...' for gth only once per protein
my $CPU = 1;                # number of CPUs that can be used
my $dir;                    # working superdirectory where programme is called from
my $genome_file;            # genome file
my $help;                   # print usage 
my $list_file;              # blast hit list file
my $log;                    # log file name
my $maxintronlen = 30000;   # maximal intron length (for usage with exonerate)
my $offset = 10000;         # offset for cutout [start-offset, end+offset] 
my $pos_file;               # position list file
my $prgsrc;                 # source programme (exonerate, spaln or gth)
my $prot_file;              # protein database file
my $prot_file_base;         # protein database file name base for $protwhole option
my %protIDs;                # hash for protein IDs
my $protWhole = 0;          # use whole protein database file
my $reg = 0;                # use regions
my %seq;                    # hash for genome sequences
my $stdoutfile;             # standard output file name
my $tmpDir;                 # temporary directory for storing protein and genome part files

my $spalnErrAdj = 80;       # version spaln2.2.0 misses a line while "counting coordinates" (+80 because of how the fasta files are printed here, see line 529-534) 
# gth options           # and other values that have proven to work well for Drosophila on chromosome 2L
my $gcmincoverage = 80; # 70 - 95 
my $prhdist  = 2;       # 2 - 6 
my $prseedlength = 20;  # 19, 20, 21 (gth: error: -prminmatchlen must be >= -prseedlength)
my $prminmatchlen = 20;

if(@ARGV==0){
  print "$usage\n"; 
  exit(0);
}

GetOptions( 'CPU=i'         => \$CPU,
            'dir=s'         => \$dir,
            'genome=s'      => \$genome_file,
            'list=s'        => \$list_file,
            'log=s'         => \$log,
            'maxintron=i'   => \$maxintronlen,
            'reg!'          => \$reg,
            'offset=i'      => \$offset,
            'pos=s'         => \$pos_file,
            'prg=s'         => \$prgsrc,
            'prot=s'        => \$prot_file,
            'help!'         => \$help);

if($help){
  print $usage;
  exit(0);
}

# if no working directory is set, use current directory
if(!defined($dir)){
 $dir = cwd(); 
}





# check whether position and list files are specified  
if(defined($pos_file)  && defined($list_file)){
  # check whether position file exists
  if(! -e $pos_file){
    print STDERR "ERROR: Position file $pos_file does not exist. Please check.\n";
    exit(1);
  }else{
    $pos_file = rel2abs($pos_file);
  }

  # check whether list file exists
  if(! -e $list_file){
    print STDERR "ERROR: List file $list_file does not exist. Please check.\n";
    exit(1);
  }else{
    $list_file = rel2abs($list_file);
  }
}else{
  $protWhole = 1;
}


# check whether protein file is specified
if(defined($prot_file)){
  # check whether protein file exists
  if(! -e $prot_file){
    print STDERR "ERROR: Protein file $prot_file does not exist. Please check.\n";
    exit(1);
  }else{
    my @a = split(/\//, $prot_file);
    $prot_file_base = $a[-1];
    $prot_file = rel2abs($prot_file);
  }
}else{
  print STDERR "ERROR: No protein file specified. Please check.\n";
  exit(1);  
}


# check whether genome file is specified
if(defined($genome_file)){
  # check whether genome file exists
  if(! -e $genome_file){
    print STDERR "ERROR: Genome file $genome_file does not exist. Please check.\n";
    exit(1);
  }else{
    $genome_file = rel2abs($genome_file);
  }
}else{
  print STDERR "ERROR: No genome file specified. Please check.\n";
  exit(1);
}



if(!defined($prgsrc)){
  print STDERR "ERROR: Please assign the source programme with --prg. Possible Options are 'exonerate', 'spaln' or 'gth'.\n";
  exit(1);
}
 
# check programme source option
if($prgsrc ne "exonerate" && $prgsrc ne "spaln" && $prgsrc ne "gth"){
  print STDERR "ERROR: Invalid value '$prgsrc' for option --prg. Possible Options are 'exonerate', 'spaln' or 'gth'. Please check.\n";
  exit(1);
}

if(!defined($log)){
  $log = "$dir/startAlign$prgsrc.log"; 
}
open (LOG, ">>".$log) or die "Cannot open file: $log\n";

$tmpDir = "$dir/tmp$prgsrc"; 
$alignDir = "$dir/align$prgsrc";

my $last_char = substr($dir, -1);
if($last_char eq "\/"){
   chop($dir);
}

if($prgsrc eq "spaln"){
# check for spaln2 environment variables
  if(!$ENV{'ALN_DBS'}){
    print STDERR "ERROR: The environment variable ALN_DBS for spaln2 is not defined. Please export an environment variable with:' export ALN_DBS=/path/to/spaln2/seqdb'\n";
    exit(1);
  }

  if(!$ENV{'ALN_TAB'}){
    print STDERR "ERROR: The environment variable ALN_TAB for spaln2 is not defined. Please export an environment variable with:' export ALN_TAB=/path/to/spaln2/table'\n";
    exit(1);
  }
}

# add 80 bp to coordinates, if --prg is spaln, 0 otherwise
if($prgsrc ne "spaln"){
  $spalnErrAdj = 0;
}


# read in list and 
if(!$protWhole){
  read_files();
}
# get protein sequences from files and store in hash
prots();

if($CPU > 1){
  get_seqs();
}
start_align();
clean_up();

close(LOG) or die("Could not close file $log!\n");


                           ############### sub functions ##############

# get protein sequences from files
sub read_files{
  open (POS, $pos_file) or die "Cannot open file: $pos_file\n";
  print LOG "\# ".(localtime).": read in positions from $pos_file\n";
  while(<POS>){
    chomp;
    my @line = split(/[\s, ]/, $_);
    if(scalar(@line) == 6){
      $contigIDs{$line[0]}{"start"} = $line[2] - 1;
      $contigIDs{$line[0]}{"end"} = $line[3] - 1;
      my @seqname = split(/[\s, ]/, $line[5]); # seqname only up to first white space
      $contigIDs{$line[0]}{"seq"} = $seqname[0];
    }
  }  
  close(POS) or die("Could not close file $pos_file!\n");

  open (LIST, $list_file) or die "Cannot open file: $list_file\n";
  print LOG "\# ".(localtime).": read in positions from $list_file\n";
  while(<LIST>){
    chomp;
    my @line = split(/[\s+, ]/, $_);
    if(scalar(@line) == 2){
      if($reg){
        push(@{$protIDs{$line[1]}}, $line[0]);
      }else{
        if(defined($contigIDs{$line[0]}{"seq"})){
          push(@{$protIDs{$line[1]}}, $contigIDs{$line[0]}{"seq"}) unless grep{$_ eq $contigIDs{$line[0]}{"seq"}} @{$protIDs{$line[1]}};
        }
      }
    }
  }  
  close(LIST) or die("Could not close file $list_file!\n");
}



sub get_seqs{
  open (FASTA, $genome_file) or die("Cannot open file: $genome_file\n");
  print LOG "\# ".(localtime).": read in DNA sequence from $genome_file\n";
  $/ = ">";
  while(<FASTA>){
    s/>$//;                           # see getAnnoFasta.pl
    next unless m/\S+/;               # see getAnnoFasta.pl
    /(.*)\n/;                         # see getAnnoFasta.pl
    my $seqname = $1;                 # see getAnnoFasta.pl
    my $sequencepart = $';            # see getAnnoFasta.pl
    $seqname =~ s/\s.*//;             # seqname only up to first white space (see getAnnoFasta.pl)
    $sequencepart =~ s/\s//g;         # see getAnnoFasta.pl
    $seq{$seqname} = $sequencepart;
  }
  $/ = "\n";
  close(FASTA) or die("Could not close file $genome_file!\n");
}


sub prots{
  if(! -d $tmpDir){
    make_path($tmpDir);
    print LOG "\# ".(localtime).": create working directory $tmpDir\n";
    print LOG "mkdir $tmpDir\n\n";
  }else{
    $cmdString = "rm -r $tmpDir";
    print LOG "\# ".(localtime).": delete existing files from $tmpDir\n";
    print LOG "$cmdString\n\n";
    system("$cmdString")==0 or die("failed to execute: $!\n");
    make_path($tmpDir);
  }

  open (PROT, $prot_file) or die "Cannot open file: $prot_file\n";
  print LOG "\# ".(localtime).": read in fasta sequences from $prot_file\n";
  my $seqname;
  my $sequencepart = "";
  while(<PROT>){
    if($_ =~ m/^>/){
      chomp;
      if($seqname){
        if($protWhole){
          $counterW++;
          print_prot("$tmpDir/$prot_file_base.addstop", $seqname, $sequencepart);
        }else{
          if(defined($protIDs{$seqname})){
            for(my $i=0; $i<scalar(@{$protIDs{$seqname}}); $i++){
              if($reg){
                if(defined($contigIDs{@{$protIDs{$seqname}}[$i]})){
                  $counterW++;
                  print_prot("$tmpDir/@{$protIDs{$seqname}}[$i].fa", $seqname, $sequencepart);
                }
              }else{
                $counterW++;
                print_prot("$tmpDir/prot@{$protIDs{$seqname}}[$i].fa", $seqname, $sequencepart);
              }
            }
          }
        }
      }
      $seqname = substr($_,1);
      $sequencepart = "";
    }else{ 
      $sequencepart .= $_;
    } 
  }

  if(defined($protIDs{$seqname})){
    for(my $i=0; $i<scalar(@{$protIDs{$seqname}}); $i++){
      if($reg){
        if(defined($contigIDs{@{$protIDs{$seqname}}[$i]})){
          print_prot("$tmpDir/@{$protIDs{$seqname}}[$i].fa", $seqname, $sequencepart);
        }
      }else{
        print_prot("$tmpDir/prot@{$protIDs{$seqname}}[$i].fa", $seqname, $sequencepart);
      }
    }
  }
  close(PROT) or die("Could not close file $prot_file!\n");
}


sub start_align{
  if(! -d $alignDir){
    make_path($alignDir);
    print LOG "\# ".(localtime).": create working directory $alignDir\n";
    print LOG "mkdir $alignDir\n\n";
  }

  chdir $tmpDir or die("Could not change to directory $tmpDir.\n");
  
  my $pm = new Parallel::ForkManager($CPU);
  my $whole_prediction_file = "$alignDir/$prgsrc.concat.aln";
  if($reg){
    foreach my $ID (keys %contigIDs){
      my $pid = $pm->start and next;
      my $target = "$contigIDs{$ID}{\"seq\"}$ID.fa";     # genome sequence
      my $query = "$ID.fa";                              # protein file
      $errorfile = "$alignDir/$contigIDs{$ID}{\"seq\"}.$ID.$prgsrc.stderr";
      $stdoutfile = "$alignDir/$contigIDs{$ID}{\"seq\"}.$ID.$prgsrc.aln";
      my $stdAdjusted;
      # create target file (contigIDs)
      if(! -e $target){ 
        my $length = $contigIDs{$ID}{"end"} - $contigIDs{$ID}{"start"} + 1 + (2 * $offset);
        my $substart = $contigIDs{$ID}{"start"} - $offset;
        my $subseq = substr($seq{$contigIDs{$ID}{"seq"}}, $substart, $length);
        print_seq($target, $subseq, $contigIDs{$ID}{"seq"});
      }
      
      if($prgsrc eq "exonerate"){
        call_exonerate($target, $query);
        $stdAdjusted = adjust_exonerate($stdoutfile, $ID);
      }
      
      if($prgsrc eq "spaln"){
        call_spaln($target, $query);
        $stdAdjusted = adjust($stdoutfile, $ID);
      }
    
      if($prgsrc eq "gth"){
        call_gth($target, $query);
        $stdAdjusted = adjust($stdoutfile, $ID);
      }
      $cmdString = "cat $stdAdjusted >>$whole_prediction_file";
      print LOG "\# ".(localtime).": add prediction from file $stdAdjusted to file $whole_prediction_file\n";
      print LOG "$cmdString\n\n";
      system("$cmdString")==0 or die("failed to execute: $!\n");
      $pm->finish;
    }
  }else{
    if($CPU > 1 || ($CPU == 1 && !$protWhole) || ($CPU == 1 && $protWhole && $prgsrc ne "gth")){
      foreach my $ID (keys %seq){
        my $pid = $pm->start and next;
        my $target = "genome$ID.fa";  # genome sequence
        my $query;
        if(!$protWhole){
          $query = "prot$ID.fa";     # protein file
        }else{
          $query = "$prot_file_base.addstop";
        }
        $errorfile = "$alignDir/$ID.$prgsrc.stderr";
        $stdoutfile = "$alignDir/$ID.$prgsrc.aln";

        # create target file (whole sequences)
        if(! -e $target){
          print_seq($target, $seq{$ID}, $ID);
        }
        if($prgsrc eq "exonerate"){
          call_exonerate($target, $query);
        }
        
        if($prgsrc eq "spaln"){
          call_spaln($target, $query);
        }
      
        if($prgsrc eq "gth"){
          call_gth($target, $query);
        }
        $cmdString = "cat $stdoutfile >>$whole_prediction_file";
        print LOG "\# ".(localtime).": add prediction from file $stdoutfile to file $whole_prediction_file\n";
        print LOG "$cmdString\n\n";
        system("$cmdString")==0 or die("failed to execute: $!\n");
        $pm->finish;
      }
    }elsif($CPU == 1 && $prgsrc eq "gth" && $protWhole){
      call_gth($genome_file, "$prot_file_base.addstop");
    }     
  }
  $pm->wait_all_children;
  chdir $dir or die ("Could not change to directory $dir.\n");
  #for version spaln2.2.0 and versions with same error
  if($prgsrc eq "spaln"){
    adjust_spaln_noreg($whole_prediction_file);
  }
}

# call assigned alignment programme
sub call_exonerate{
  my $tFile = shift;
  my $qFile = shift;
  $cmdString = "exonerate --model protein2genome --maxintron $maxintronlen --showtargetgff --showalignment no --query $qFile -t   $tFile > $stdoutfile 2>$errorfile";
  print LOG "\# ".(localtime).": run exonerate for target '$tFile' and query '$qFile'\n";
  print LOG "$cmdString\n\n";
  system("$cmdString")==0 or die("failed to execute: $!\n");
}

sub call_gth{
  my $genomic = shift;
  my $protein = shift;
  $cmdString = "gth -genomic $genomic -protein $protein -gff3out -skipalignmentout -paralogs -prseedlength $prseedlength -prhdist $prhdist -gcmincoverage $gcmincoverage -prminmatchlen $prminmatchlen -o  $stdoutfile 2>$errorfile";
  print LOG "\# ".(localtime).": run GenomeThreader for genome '$genomic' and protein '$protein'\n";
  print LOG "$cmdString\n\n";
  system("$cmdString")==0 or die("failed to execute: $!\n");
}

sub call_spaln{
  my $tFile = shift;
  my $qFile = shift;
  my @split = split(/\./, $tFile);
  if(! -f "$split[0].idx"){ 
    $cmdString = "makdbs -KP $tFile;";
    print LOG "\# ".(localtime).": create *.ent, *.grp, *.idx, (*.odr), and *.seq files for target '$tFile' \n";
    print LOG "$cmdString\n\n";
    system("$cmdString")==0 or die("failed to execute: $!\n");
  }
  if(! -f "$split[0].bkp"){
    $cmdString = "perl $ENV{'ALN_DBS'}/makblk.pl -W$split[0].bkp -KP $tFile";
    print LOG "\# ".(localtime).": create *.bkp file for target '$tFile'\n";
    print LOG "$cmdString\n\n";
    system("$cmdString")==0 or die("failed to execute: $!\n");
  }

  @split = split(/\./, $qFile);
  if(! -f "$split[0].idx"){
    $cmdString = "makdbs -KA $qFile;";
    print LOG "\# ".(localtime).": create *.ent, *.grp, *.idx, (*.odr), and *.seq files for query '$qFile' \n";
    print LOG "$cmdString\n\n";
    system("$cmdString")==0 or die("failed to execute: $!\n");
  }
  if(! -f "$split[0].bka"){
    $cmdString = "perl $ENV{'ALN_DBS'}/makblk.pl -W$split[0].bka -KA $qFile";
    print LOG "\# ".(localtime).": create *.bka file for query '$qFile'\n";
    print LOG "$cmdString\n\n";
    system("$cmdString")==0 or die("failed to execute: $!\n");
  }

  $cmdString = "spaln -Q7 -O0 $tFile $qFile > $stdoutfile 2>$errorfile";
  print LOG "\# ".(localtime).": run spaln for target '$tFile' and query '$qFile'\n";
  print LOG "$cmdString\n\n";
  system("$cmdString")==0 or die("failed to execute: $!\n");
}

sub print_prot{
  my $outfile = shift;
  my $Sname = shift;
  my $seqpart = shift;
  open(OUT, ">>$outfile") or die "Cannot open file: $outfile\n";
  print OUT ">$Sname\n";
  # add '*' for GenomeThreader (gth) [instead of 'gt seqtransform -addstopaminos ...' since this is not part 
  # of gth, but of GenomeTools (gt)] 
  if($prgsrc eq "gth"){
    print OUT "$seqpart*\n";
    if($counterW == 1){
      print LOG "\# ".(localtime).": add stop amino acid ('*') to protein sequences'\n";
    }
  }else{
    print OUT "$seqpart\n";
  }
  close(OUT) or die("Could not close file $outfile!\n");
}


# print genome sequence part or whole sequence
sub print_seq{
  my $tFile = shift;
  my $sequence = shift;
  my $Sname = shift;
  open(TARGET, ">$tFile") or die "Cannot open file: $tFile\n";
  print TARGET ">$Sname\n";
  my $start = 0; # see getAnnoFasta.pl
  my $ret = "";  # see getAnnoFasta.pl
  while (length($sequence)-$start >= 80) {     # see getAnnoFasta.pl
    my $shortline = substr($sequence, $start, 80);
    $ret .= "$shortline\n";
    $start += 80;
  }
  $ret .= substr($sequence, $start, 80) . "\n" if ($start<length($sequence));
  print TARGET $ret;
  close(TARGET) or die("Could not close file $tFile!\n");
}


# adjust spaln and gth gff3 output
sub adjust{
  my $output = shift;
  my $conID = shift;
  my $output_adj = "$output.adj";
  open(OUT, ">$output_adj") or die "Cannot open file: $output_adj\n";
  open(IN, $output) or die "Cannot open file: $output\n";
  while(<IN>){
    chomp;
    my @f = split(/\s/, $_);
    if($_ =~ m/##sequence-region/){
        $f[-1] += $contigIDs{$conID}{"start"} - $offset + $spalnErrAdj;
        $f[-2] += $contigIDs{$conID}{"start"} - $offset + $spalnErrAdj; 
      print OUT "##sequence-region\t$f[-3] $f[-2] $f[-1]\n";
    }elsif(scalar(@f) >= 9){
      $f[3] += $contigIDs{$conID}{"start"} - $offset + $spalnErrAdj;
      $f[4] += $contigIDs{$conID}{"start"} - $offset + $spalnErrAdj;
      print OUT "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\t$f[8]\n";
    }else{
      print OUT "$_\n";
    }
  }
  close(IN) or die("Could not close file $output!\n");
  close(OUT) or die("Could not close file $output_adj!\n");
  return $output_adj;
}

# adjust exonerate output
sub adjust_exonerate{
  my $exonerate_output = shift;
  my $conID = shift;
  my $exonerate_output_adj = "$exonerate_output.adj";
  open(OUT, ">$exonerate_output_adj") or die "Cannot open file: $exonerate_output_adj\n";
  open(IN, $exonerate_output) or die "Cannot open file: $exonerate_output\n";
  my $add_space = "";
  while(<IN>){
    chomp;
    if($_ =~ m/  Target range: (\d+) -> (\d+)/){
      my $start = $1 + $contigIDs{$conID}{"start"} - $offset;
      my $end = $2 + $contigIDs{$conID}{"start"} - $offset;
      my $diff = length($end) - length($2);
      $add_space = " " x $diff; 
      print OUT "  Target range: $start -> $end\n";
    }elsif($_ =~m/^(\s+)(\d+) : ([^a-z\.\-]+) : (\d+)/){
      my $start = $2 + $contigIDs{$conID}{"start"} - $offset;
      my $end = $4 + $contigIDs{$conID}{"start"} - $offset;
      print OUT "$1"."$start : $3 : "."$end\n";
    }elsif($_ =~ m/vulgar:/){
      my @line = split(/[\s, ]/, $_);
      $line[6] += $contigIDs{$conID}{"start"} - $offset;
      $line[7] += $contigIDs{$conID}{"start"} - $offset;
      my $vulgar = join(' ', @line);
      print OUT "$vulgar\n";
    }elsif($_ =~m/^(\s+)(\d+) : ([a-zA-Z\.\-]+)/ || $_ =~m/^(\s+)([a-zA-Z\.\-]+)/ || $_ =~ m/^(\s+)[\|\.\!\:]+/){
      print OUT $add_space.$_."\n";
    }else{
      my @f = split(/\t/, $_);
      if(scalar(@f) == 9 || scalar(@f) == 8){
        if($f[6] eq "+"){
          $f[3] += $contigIDs{$conID}{"start"} - $offset;
          $f[4] += $contigIDs{$conID}{"start"} - $offset;
        }else{
          $f[3] += $contigIDs{$conID}{"start"} - $offset;
          $f[4] += $contigIDs{$conID}{"start"} - $offset;
        }
        print OUT "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]";
        if(scalar(@f) == 9){  
          if($f[2] eq "similarity"){
            my @sim_f = split(/\t/, $f[8]);
            for(my $i = 0; $i < scalar(@sim_f); $i++){
              if($sim_f[$i] eq "Align"){
                $sim_f[$i+1] += $contigIDs{$conID}{"start"};
              }
            }
            $f[8] = join(' ', @sim_f);
          }
          print OUT "\t$f[8]";
        }
        print OUT "\n";
      }else{
        print OUT "$_\n";
      }
    }
  }
  close(IN) or die("Could not close file $exonerate_output!\n");
  close(OUT) or die("Could not close file $exonerate_output_adj!\n");
  return $exonerate_output_adj;
}


# adjust spaln output for noreg option
sub adjust_spaln_noreg{
  my $output = shift;
  my $output_adj = "$output.adj";
  open(OUT, ">$output_adj") or die "Cannot open file: $output_adj\n";
  open(IN, $output) or die "Cannot open file: $output\n";
  while(<IN>){
    chomp;
    my @f = split(/\s/, $_);
    if($_ =~ m/##sequence-region/){
        $f[-1] += $spalnErrAdj;
        $f[-2] += $spalnErrAdj; 
      print OUT "##sequence-region\t$f[-3] $f[-2] $f[-1]\n";
    }elsif(scalar(@f) >= 9){
      $f[3] += $spalnErrAdj;
      $f[4] += $spalnErrAdj;
      print OUT "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\t$f[8]\n";
    }else{
      print OUT "$_\n";
    }
  }
  close(IN) or die("Could not close file $output!\n");
  close(OUT) or die("Could not close file $output_adj!\n");
}

# delete empty files
sub clean_up{
  print STDOUT "NEXT STEP: delete empty files\n";
  my @files = `find $alignDir -empty`;
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
