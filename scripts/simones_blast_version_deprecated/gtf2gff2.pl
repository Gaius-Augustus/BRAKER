#!/usr/bin/perl
#
# format convert a gtf file
#
# Mario Stanke, 1.2.2010, mario.stanke@uni-greifswald.de
use strict;
use Getopt::Long;


my $help = 0;
my $verbose = 0;
my $printExon = 0;
my $printUTR = 0;
my $gff3 = 0;
my $printIntron = 0;
my $outfile;
my $genemarkET = 0; # added by Simone Lange: only necessary for GeneMark-ET versions without exon entries because otherwise one will get this error: Can't use an undefined value as an ARRAY reference at gtf2gff.pl line 217, <STDIN> line 292186. with --printIntron or Can't use an undefined value as an ARRAY reference at gtf2gff.pl line 265, <STDIN> line 292186. without it




GetOptions('out=s'=>\$outfile,
           'help!'=>\$help,
           'verbose!'=>\$verbose,
           'printExon!'=>\$printExon,
           'printUTR!'=>\$printUTR,
           'printIntron!'=>\$printIntron,
           'gff3!'=>\$gff3,
           'GMET!'=>\$genemarkET);

exec("perldoc $0") if ($help || !defined($outfile));

my @txorder = (); # tx ids in the order of the input file
my %geneOf = ();  # keys tx ids, values: gene ids
my %geneLine =  (); # keys gene ids, values: array refs for the gene GTF line (if exists)
# hash of transcripts
#    keys: transcript ids
#    values: hash reference
#            keys: txstart
#                  txend
#                  codingstart
#                  codingend
#                  strand, chr, source
#                  txline  array of columns if a transcript/mRNA line exists
#                  CDS  array of arrays (lines and columns) for coding parts of exons
#                  UTR  array of arrays for UTR exons
#                  exon array of arrays for complete exons
#                  intron array of arrays for introns
#                  rest array of arrays all other features, like tts,tss start_codon, stop_codon

my %txs = ();

parseAndStoreGTF();
convert();
open (OUT, ">$outfile") or die ("Could not open $outfile for writing.");
printConvertedGTF();
close OUT;



sub parseAndStoreGTF{
  my %seen = ();
  my ($txid, $geneid, $chr, $start, $end, $feature, $strand, $source);
  foreach my $line (<STDIN>){
    my @f = split /\t/, $line;
    next if (@f<8);
    ($chr,$source,$feature,$start,$end,$strand) = ($f[0],$f[1],$f[2],$f[3],$f[4],$f[6]);
    # check whether it is a line with 'gene' feature
    if ($f[2] eq "gene" && ($f[8] =~ /^(\S+)$/ || $f[8] =~ /gene_id."?([^";]+)"?/)){
      @_=split(/\./,$1);        # changed by Simone Lange for use with wormbase.gtf
      $geneid =join("_",@_);    # changed by Simone Lange for use with wormbase.gtf, previously $geneid = $1;
      $geneLine{$geneid} = \@f;
      next;
    } 
    # check whether it is a line with 'transcript' feature
    if ($f[2] eq "transcript" && ($f[8] =~ /^(\S+)$/ || $f[8] =~ /transcript_id."?([^";]+)"?/)){
      @_=split(/\./,$1);        # see above
      $txid =join("_",@_).".1"; # see above, previously $txid = $1;
      $txs{$txid} = {"strand"=>$strand, "chr"=>$chr, "source"=>$source, "CDS"=>[], "UTR"=>[], "exon"=>[], "intron"=>[], "rest"=>[]} if (!exists($txs{$txid}));
      $txs{$txid}{"txline"} = \@f;
      next;
    }
        
    $txs{$txid}{"CDS"} = [] if (!defined($txs{$txid}{"CDS"}));
    $txs{$txid}{"UTR"} = [] if (!defined($txs{$txid}{"UTR"}));
    $txs{$txid}{"exon"} = [] if (!defined($txs{$txid}{"exon"}));
    $txs{$txid}{"rest"} = [] if (!defined($txs{$txid}{"rest"}));

    # all other lines must belong to a transcript and a gene
    if ($f[8] =~ /(transcript_id|Transcript)."?([^";]+)"?/){
      @_=split(/\./,$2);         # see above
      $txid = join("_",@_).".1"; # see above, previously $txid = $2;
    } else {
      die ("Not GTF format in the following line:\n$line\ntranscript_id not found.\n");
    }
    if ($f[8] =~ /gene_id."?([^";]+)"?/){
      @_=split(/\./,$1);     # see above
      $geneid =join("_",@_); # see above, previously $geneid = $1;
    } else {
      die ("Not GTF format in the following line:\n$line\ngene_id not found.\n");
    }
    if (!$seen{$txid}){
      push @txorder, $txid; # remember the input order for transcripts for the output
      $seen{$txid} = 1;
    }
    # assign parental gene id to tx id
    die ("transcript $txid has conflicting gene parents: and $geneid. Remember: In GTF txids need to be overall unique.")
    if (defined($geneOf{$txid}) && $geneOf{$txid} ne $geneid);
      $geneOf{$txid} = $geneid;

    if ($feature eq "CDS" || $feature eq "coding_exon" || $feature eq "exon" || $feature =~ /UTR/i){
      $txs{$txid} = {"strand"=>$strand, "chr"=>$chr, "source"=>$source, "CDS"=>[], "UTR"=>[], "exon"=>[], "intron"=>[], "rest"=>[]} if (!exists($txs{$txid}));
      $txs{$txid}{"txstart"} = $start if (!defined($txs{$txid}{"txstart"}) || $txs{$txid}{"txstart"} > $start);
      $txs{$txid}{"txend"} = $end if (!defined($txs{$txid}{"txend"}) || $txs{$txid}{"txend"} < $end);
    }
    if ($feature eq "CDS" || $feature eq "coding_exon"){
      $txs{$txid}{"codingstart"} = $start if (!defined($txs{$txid}{"codingstart"}) || $txs{$txid}{"codingstart"} > $start);
      $txs{$txid}{"codingend"} = $end if (!defined($txs{$txid}{"codingend"}) || $txs{$txid}{"codingend"} < $end);
      push @{$txs{$txid}{"CDS"}}, \@f;
    } elsif ($feature =~ /UTR/i){
      push @{$txs{$txid}{"UTR"}}, \@f;
    } elsif ($feature eq "exon"){
      push @{$txs{$txid}{"exon"}}, \@f;
    } elsif ($feature eq "intron") {
      $txs{$txid}{"intron"} = [] if (!defined($txs{$txid}{"intron"}));
      push @{$txs{$txid}{"intron"}}, \@f;
    } else {
      push @{$txs{$txid}{"rest"}}, \@f;
    }
  }
}


sub convert{
  my @f;
  foreach my $txid (keys %txs){
  # remember whether exons were not in the input file
    my $exonArrayWasEmpty;
    if(@{$txs{$txid}{"exon"}} == 0){
      $exonArrayWasEmpty = 1;
    }else{
      $exonArrayWasEmpty = 0;
    }
    # add exon lines if not already present and if desired in output
    if (($printExon || $printIntron) && @{$txs{$txid}{"exon"}} == 0){
      print "Creating exon lines for $txid\n" if ($verbose);
      # sort UTR and CDS lines by coordinates
      my @exonpartlines = sort {$a->[3] <=> $b->[3] || $a->[4] <=> $b->[4]} (@{$txs{$txid}{"CDS"}}, @{$txs{$txid}{"UTR"}});
      next if (@exonpartlines == 0);
      @f = @{$exonpartlines[0]};
      shift @exonpartlines;
      ($f[2], $f[5], $f[7]) = ("exon", '.', '.'); # score and frame are not defined
      foreach my $g (@exonpartlines){
        if ($f[4] >= $g->[3]){ # check for non-overlappingness
          die ("In transcript $txid two UTR/CDS features are overlapping. Not allowed by definition.");
        } elsif ($f[4] + 1 == $g->[3]){ # exactly adjacent
          # join two UTR/CDS features to one
          $f[4] = $g->[4];
        } else {
          # push exon
          my @ff = @f; # deep copy array
          push @{$txs{$txid}{"exon"}}, \@ff;
          @f = @$g;
          ($f[2], $f[5], $f[7]) = ("exon", '.', '.'); # score and frame are not defined
        }
      }
      # push remaining, last exon
      my @ff = @f;
      push @{$txs{$txid}{"exon"}}, \@ff;
    }

    # add UTR lines if not already present and if desired in output
    if ($printUTR && @{$txs{$txid}{"UTR"}} == 0){
      # sort exon and CDS lines start coordinates, "exon" feature comes before CDS feature at tie
      my @epl = sort { return 1 if ($a->[3] > $b->[3]);
                       return -1 if ($a->[3] < $b->[3]);
                       return 1 if ($a->[2] ne "exon");
                       return -1 if ($b->[2] ne "exon");
                       return 0;} (@{$txs{$txid}{"CDS"}}, @{$txs{$txid}{"exon"}});
      for (my $i=0; $i<@epl; $i++){
        next if ($epl[$i]->[2] ne "exon");
        if ($i < @epl-1 && $epl[$i+1]->[2] ne "exon" && $epl[$i+1]->[3] <= $epl[$i]->[4]){ # next feature is CDS and overlapping somehow
          # --------------- exon
          #    XXXXXXX CDS
          # left UTR part
          if ($epl[$i]->[3] < $epl[$i+1]->[3]){
            my @ff = @{$epl[$i]};
            $ff[4] = $epl[$i+1]->[3]-1; # ends 1 before CDS starts
            ($ff[5], $ff[7]) = ('.', '.'); # score and frame are not defined
            $ff[2] = ($txs{$txid}{"strand"} eq '+')? "5'-UTR" : "3'-UTR";
            push @{$txs{$txid}{"UTR"}}, \@ff;
          }
          # right UTR part
          if ($epl[$i+1]->[4] < $epl[$i]->[4]){
            my @ff = @{$epl[$i]};
            $ff[3] = $epl[$i+1]->[4]+1; # starts 1 after CDS starts
            ($ff[5], $ff[7]) = ('.', '.'); # score and frame are not defined
            $ff[2] = ($txs{$txid}{"strand"} eq '+')? "3'-UTR" : "5'-UTR";
            push @{$txs{$txid}{"UTR"}}, \@ff;
          }
        } else {
          my @ff = @{$epl[$i]};
          ($ff[5], $ff[7]) = ('.', '.'); # score and frame are not defined
          if (($ff[3] < $txs{$txid}{"codingstart"} && $txs{$txid}{"strand"} eq "+") ||
            ($ff[4] > $txs{$txid}{"codingend"} && $txs{$txid}{"strand"} eq "-")){
            $ff[2] = "5'-UTR";
          } else {
            $ff[2] = "3'-UTR";
          }
          push @{$txs{$txid}{"UTR"}}, \@ff;
        }
      }
    }
    # add intron lines if not already present and if desired in output 
    # Simone Lange: adjusted if-condition for compatibility with GeneMark-ET input; previously  
    #if ($printIntron && @{$txs{$txid}{"intron"}} == 0){     
    if (($printIntron && !defined(@{$txs{$txid}{"intron"}})) || ($printIntron && @{$txs{$txid}{"intron"}} == 0)){
      my @exonlines =  sort {$a->[3] <=> $b->[3] || $a->[4] <=> $b->[4]} (@{$txs{$txid}{"exon"}});
      my @intronline = @{$exonlines[0]};
      my $start = $intronline[4]+1;
      my $end;
      shift @exonlines;
      foreach my $ex (@exonlines){
        $end = $ex->[3]-1;
        my @il = @intronline;
        ($il[2], $il[3], $il[4], $il[5], $il[7]) = ("intron",$start,$end,'.','.');
        push @{$txs{$txid}{"intron"}}, \@il;
        $start = $ex->[4]+1;
      }   
      if(!$printExon && $exonArrayWasEmpty){
        @{$txs{$txid}{"exon"}} = ();
                
      }
    }
  }
}


sub printConvertedGTF {
  my @lines;
  my %seen = ();
  my $geneid;
  foreach my $txid (@txorder){
    # print gene line before the first transcript of this gene
    $geneid = $geneOf{$txid};
    if (!$seen{$geneid} && defined($geneLine{$geneid})){
      if ($gff3) {
        $geneLine{$geneid}->[8] = "ID=$geneid;\n";
      }
      print OUT join ("\t", @{$geneLine{$geneid}});
      $seen{$geneid} = 1;
    }
    # print transcript line
    if ($txs{$txid}{"txline"}[0] ne ""){
      if ($gff3) {
        $txs{$txid}{"txline"}->[2] = "mRNA";
        $txs{$txid}{"txline"}->[8] = "ID=$txid;Parent=$geneid\n";
      }
      print OUT join ("\t", @{$txs{$txid}{"txline"}});
    }
    # print all other lines
    # Simone Lange: remove @{$txs{$txid}{"intron"}}, when genemark.gtf does not contain exon entries and --printIntron option is not used; previously without if-condition 
    if($genemarkET && !$printIntron){ 
      @lines = sort {$a->[3] <=> $b->[3] || $a->[4] <=> $b->[4]}
      (@{$txs{$txid}{"CDS"}}, @{$txs{$txid}{"UTR"}}, @{$txs{$txid}{"exon"}}, @{$txs{$txid}{"rest"}}); 
    }else{
      @lines = sort {$a->[3] <=> $b->[3] || $a->[4] <=> $b->[4]}
      (@{$txs{$txid}{"CDS"}}, @{$txs{$txid}{"UTR"}}, @{$txs{$txid}{"exon"}}, @{$txs{$txid}{"intron"}}, @{$txs{$txid}{"rest"}});
    }
        #[PKRK] following variable are to make the CDS and exon ID unique within the scope of the gene model
    my $ct_exon = 0;
    my $ct_CDS = 0;
    my $ct_3UTR = 0;
    my $ct_5UTR = 0;
    foreach my $line (@lines){
      if ($gff3){
        $line->[2] =~ s/3.*UTR/three_prime_utr/;
        $line->[2] =~ s/5.*UTR/five_prime_utr/;
        $line->[2] =~ s/tss/transcription_start_site/;
        $line->[2] =~ s/tts/transcription_end_site/;
        if($line->[2] eq "exon"){
          ++$ct_exon;
          $line->[8] = "ID=$txid.$line->[2]$ct_exon;Parent=$txid;\n";
        }elsif($line->[2] eq "CDS"){
          ++$ct_CDS;
          $line->[8] = "ID=$txid.$line->[2]$ct_CDS;Parent=$txid\n";
        }elsif($line->[2] eq "three_prime_utr"){
          ++$ct_3UTR;
          $line->[8] = "ID=$txid.3UTR$ct_3UTR;Parent=$txid\n";
        }elsif($line->[2] eq "five_prime_utr"){
          ++$ct_5UTR;
          $line->[8] = "ID=$txid.5UTR$ct_5UTR;Parent=$txid\n";
                
        } else {
          $line->[8] = "Parent=$txid;\n";
        }
      }
      print OUT join("\t", @$line);
    }
  }
}


__END__

=pod

=head1 NAME

gtf2gff.pl      format convert a GTF file

=head1 SYNOPSIS

gtf2gff.pl <in.gtf --out=out.gff

    Besides easy by-line changes this script can in particular swap between the different representations of UTRs:
    a) explicit UTR lines (e.g. 3'-UTR or three_prime_utr)
    b) implicit specification by 'exon' and 'CDS' features
    
=head1 OPTIONS

  out              gff output file
  --printExon      print exon features (may include CDS and UTR parts)
  --printUTR       print UTR features
  --printIntron    print intron features
  --gff3           output in gff3 format

=head1 DESCRIPTION
    
    example input:

    chr1 AUGUSTUS  gene        12656   14013   0.04    +   .   g50
    chr1 AUGUSTUS  transcript  12656   14013   0.04    +   .   g50.t1
    chr1 AUGUSTUS  tss         12656   12656   .       +   .   transcript_id "g50.t1"; gene_id "g50";
    chr1 AUGUSTUS  5'-UTR      12656   12867   0.2     +   .   transcript_id "g50.t1"; gene_id "g50";
    chr1 AUGUSTUS  start_codon 12868   12870   .       +   0   transcript_id "g50.t1"; gene_id "g50";
    chr1 AUGUSTUS  intron      12994   13248   1       +   .   transcript_id "g50.t1"; gene_id "g50";
    chr1 AUGUSTUS  CDS         12868   12993   0.8     +   0   transcript_id "g50.t1"; gene_id "g50";
    chr1 AUGUSTUS  CDS         13249   13479   1       +   0   transcript_id "g50.t1"; gene_id "g50";
    chr1 AUGUSTUS  stop_codon  13477   13479   .       +   0   transcript_id "g50.t1"; gene_id "g50";
    chr1 AUGUSTUS  3'-UTR      13480   14013   0.17    +   .   transcript_id "g50.t1"; gene_id "g50";
    chr1 AUGUSTUS  tts         14013   14013   .       +   .   transcript_id "g50.t1"; gene_id "g50";

    example output for --gff3 --printExon:

    chr1 AUGUSTUS  gene                     12656   14013   0.04    +   .   ID=g50;
    chr1 AUGUSTUS  mRNA                     12656   14013   0.04    +   .   ID=g50.t1;Parent=g50
    chr1 AUGUSTUS  transcription_start_site 12656   12656   .       +   .   Parent=g50.t1;
    chr1 AUGUSTUS  five_prime_utr           12656   12867   0.2     +   .   ID=g50.t1.5UTR1;Parent=g50.t1
    chr1 AUGUSTUS  exon                     12656   12993   .       +   .   ID=g50.t1.exon1;Parent=g50.t1;
    chr1 AUGUSTUS  start_codon              12868   12870   .       +   0   Parent=g50.t1;
    chr1 AUGUSTUS  CDS                      12868   12993   0.8     +   0   ID=g50.t1.CDS1;Parent=g50.t1
    chr1 AUGUSTUS  intron                   12994   13248   1       +   .   Parent=g50.t1;
    chr1 AUGUSTUS  CDS                      13249   13479   1       +   0   ID=g50.t1.CDS2;Parent=g50.t1
    chr1 AUGUSTUS  exon                     13249   14013   .       +   .   ID=g50.t1.exon2;Parent=g50.t1;
    chr1 AUGUSTUS  stop_codon               13477   13479   .       +   0   Parent=g50.t1;
    chr1 AUGUSTUS  three_prime_utr          13480   14013   0.17    +   .   ID=g50.t1.3UTR1;Parent=g50.t1
    chr1 AUGUSTUS  transcription_end_site   14013   14013   .       +   .   Parent=g50.t1;

=cut
