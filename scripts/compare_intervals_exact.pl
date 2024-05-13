#!/usr/bin/env perl
# ==============================================================
# Alex Lomsadze, Tomas Bruna
# Georgia Institute of Technology, Atlanta, Georgia, US
# Last update October 2023
# 
# This script compares intervals from GFF/GTF/GFF3-like files and
# generates a report suitable for drawing a Venn diagram.
#
# Comparison can be done between 2 or 3 annotations.
#
# Many non-default output options are applicable only to input files in enriched GeneMark GTF format
# Optional parameters can be used to account for pseudogenes, repeat regions, etc.
# ==============================================================

use strict;
use warnings;

use Getopt::Long qw( GetOptions );
use Storable qw(dclone);
use Data::Dumper;

my $VERSION = "v11_2023";

# ------------------------------------------------
my $v = 0;     # verbose
my $debug = 0;

my $f1 = '';   # input file 1 with coding genes in GFF-like format
my $f2 = '';   # input file 2 with coding genes in GFF-like format 
my $f3 = '';   # input file 3 with coding genes in GFF-like format

my $pseudo = '';   # input file with pseudogene coordinates in GFF-like format
my $masked = '';   # input file with repeat coordinates in GFF-like format
my $min_masked = 1000;   # ignor repeats shorter than minimum

my $compare_cds = 0;  # default comparision mode; value "1" is set in ParseCMD()
my $compare_introns = 0;
my $compare_donors = 0;
my $compare_acceptors = 0;
my $compare_starts = 0;
my $compare_stops = 0;
my $compare_single = 0;
my $compare_initial = 0;
my $compare_internal = 0;
my $compare_terminal = 0;
my $compare_gene = 0;
my $compare_multiGene = 0;   # multi-cds genes only
my $compare_singleGene = 0;  # single-cds genes only
my $compare_trans = 0;
my $compare_exon = 0;

my $no_phase = '';  # ignore phase of feature in comparision

my $out_file = '';
my $original = 0;

my $shared12 = '';
my $shared13 = '';
my $shared23 = '';

my $shared123 = '';

my $unique1 = '';
my $unique2 = '';
my $unique3 = '';
# ------------------------------------------------

Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();

# hashes "%h*" have keys formed from interval coordinates
# for example, GFF fields 1,4,5,7 and 8 are joined together with underscores: seqid_start_end_strand_phase

my %h1;
my %h2;
my %h3;

# hashes "%tr2gene*" are used only for gene level accuracy calculation
# transcript id - to - gene id hash
# in the case of gene with isoforms : transcripts ids (keys) will have the same value in hash - gene ID

my %tr2gene1 = ();
my %tr2gene2 = ();
my %tr2gene3 = ();

# parse GFF files with annotation

ParseGFF( $f1, \%h1, \%tr2gene1 );
ParseGFF( $f2, \%h2, \%tr2gene2 );
ParseGFF( $f3, \%h3, \%tr2gene3 ) if $f3; 

# Remove intervals from sets if they overlap with pseudo or mask
# Removal rules differ for pseudo and mask

FilterPseudo( $pseudo ) if $pseudo;
FilterMasked( $masked, $min_masked ) if $masked;

# Gene-level comparison logic differs from the logic used in the comparison of all the other interval types
# Gene level comparison is done in the first IF below
# Gene level results are then placed into %h1 and %h2 hashes
# The same interface and output style are used later for genes as for other features.

if ( $compare_gene or $compare_multiGene or $compare_singleGene )
{
	if ($f3)
		{ die "sorry, 3-way gene level comparision is not implemented yet\n"; }

	# %h* hash key -> to array of trans_id's
	# replace trans_id by gene_id
	# %h* hash key -> to gene_id
	# multy trans_id is sorted - first one is used

	ReplaceValue( \%h1, \%tr2gene1 );
	ReplaceValue( \%h2, \%tr2gene2 );

	CountGenes(\%tr2gene1, $f1) if $v;
	CountGenes(\%tr2gene2, $f2) if $v;

	if ($debug)
	{
		CountDistrTrPerGene(\%h1);
		CountDistrTrPerGene(\%h2);
	}

	# compare two gene sets
	# comparison results are saved into %z* hashes
	# %z* hashes - keys gene_id

	# put all gene_id's from %h1 into %z1
	my %z1;
	foreach my $key (keys %h1)
	{
		my $gid_new = $h1{$key} ." f1";
		$z1{ $gid_new } += 1;
	}

	# put all gene_id's from %h2 into %z2
	my %z2;
	foreach my $key (keys %h2)
	{
		my $gid_new = $h2{$key} ." f2";
		$z2{$gid_new} += 1;
	}

	# put matching genes into hash
	my %z1_z2_overlap = ();
	my %z2_z1_overlap = ();

	foreach my $key (keys %h1)
	{
		if ( exists $h2{$key} )
		{
			my $gid_new = $h1{$key} ." f1";
			$z1_z2_overlap{ $gid_new } = $h2{$key} ." f2";
		}
	}

	foreach my $key (keys %h2)
	{
		if ( exists $h1{$key} )
		{
			my $gid_new = $h2{$key} ." f2";
			$z2_z1_overlap{ $gid_new } = $h1{$key} ." f1";
		}
	}

	if ($v)
	{
		print "# in z1: ". (scalar (keys %z1)) ."\n";
		print "# in z2: ". (scalar (keys %z2)) ."\n";
		print "# in z1-z2: ". (scalar (keys %z1_z2_overlap)) ."\n";
		print "# in z2-z1: ". (scalar (keys %z2_z1_overlap)) ."\n";
	}

	# replace gene_id in %z2 by name from %z1 for %z2_z1_overlap
	foreach my $key (keys %z2_z1_overlap)
	{
		if ( exists $z2{$key} )
		{
			delete $z2{$key};
			$z2{ $z2_z1_overlap{$key} } += 1;
		}
	}

	%h1 = %z1;
	%h2 = %z2;

	if ($v)
	{
		print "# in z2: ". (scalar (keys %z2)) ."\n";
	}
}

# variables below hold information from comparisons

my %shared_12;
my %shared_13;
my %shared_23;

my %shared_123;

my %unique_1;
my %unique_2;
my %unique_3;

# run comparison

if( !$f3 )
{
	Compare2();
}
else
{
	Compare3();
}

exit 0;

# ================= subs =========================
sub CountGenes
{
	my $ref = shift;
	my $label = shift;

	my %genes = ();

	foreach my $tid (keys %{$ref})
	{
		$genes{ $ref->{$tid} } += 1;
	}

	print "# genes in set $label: ". (scalar keys %genes) ."\n"; 
}
# ------------------------------------------------
sub CountDistrTrPerGene
{
	my $ref = shift;

	my %gene_count;

	foreach my $key (keys %{$ref})
	{
		$gene_count{$ref->{$key}} += 1;
	}

	my %hist;

	foreach my $key (keys %gene_count)
	{
		$hist{ $gene_count{$key} } += 1;
	}

	print "# transc-per-gene hist\n";

	foreach my $key ( sort{$a<=>$b} keys %hist )
	{
		print "# $key $hist{$key}\n";
	}
}
# ------------------------------------------------
sub ReplaceValue
{
	my $ref = shift;
	my $new_v = shift;

	foreach my $key ( keys %{$ref} )
	{
		my $trans_id = '';

		if( $ref->{$key} =~ /^(\S+)\s*/ )
		{
			$trans_id = $1;
		}
		else
			{ die "error, unexpected format for key: $key\n"; }

		if ( exists $new_v->{$trans_id} )
		{
			$ref->{$key} = $new_v->{$trans_id};
		}
		else
			{ die "error, gene id is missing for transcript id: $trans_id\n"; }
	}
}
# ------------------------------------------------
sub FilterMasked
{
	my $name = shift;
	my $min = shift;

	my %masked_h = LoadIntervals( $name, $min );

	my %intervals = ();
	my $label = '';

	%intervals = HashKeysToArrIntervals( \%h1 );
	$label = "# intervals removed in --f1 due to overlap with masked regions: ";
	RemoveOverlapping( \%masked_h, \%intervals, 1, $label, \%h1 );

	%intervals = HashKeysToArrIntervals( \%h2 );
	$label = "# intervals removed in --f2 due to overlap with masked regions: ";
	RemoveOverlapping( \%masked_h, \%intervals, 1, $label, \%h2 );

	if ( $f3 )
	{
		%intervals = HashKeysToArrIntervals( \%h3 );
		$label = "# intervals removed in --f3 due to overlap with masked regions: ";
		RemoveOverlapping( \%masked_h, \%intervals, 1, $label, \%h3 );
	}
}
# ------------------------------------------------
sub FilterPseudo
{
	my $name = shift;

	my %pseudo_h = LoadIntervals( $name, 0 );

	my %intervals = ();
	my $label = '';

	%intervals = HashKeysToArrIntervals( \%h1 );
	$label = "# pseudo regions removed due to overlap with --f1: ";
	RemoveOverlapping( \%pseudo_h, \%intervals, 0, $label );
	
	%intervals = HashKeysToArrIntervals( \%h2 );
	$label = "# intervals removed in --f2 due to overlap with pseudo regions: ";
	RemoveOverlapping( \%pseudo_h, \%intervals, 1, $label, \%h2 );

	if ( $f3 )
	{
		%intervals = HashKeysToArrIntervals( \%h3 );
		$label = "# intervals removed in --f2 due to overlap with pseudo regions: ";
		RemoveOverlapping( \%pseudo_h, \%intervals, 1, $label, \%h3 );
	}
}
# ------------------------------------------------
sub RemoveOverlapping
{
	my $aa = shift;
	my $bb = shift;
	my $in_bb = shift;
	my $label = shift;
	my $ref = shift;

	my @arr = ();
	my %hash = ();

	foreach my $key ( keys %{$aa} )
	{
		next if ( ! exists $bb->{$key} );

		my $size_aa = scalar @{$aa->{$key}};
		my $size_bb = scalar @{$bb->{$key}};

		my $i = 0;
		my $j = 0;
		
		print "# resolving overlaps for $key $size_aa $size_bb\n" if $debug;

		while( $i < $size_aa && $j < $size_bb )
		{
			# i.R < j.L
			if ( $aa->{$key}[$i][1] < $bb->{$key}[$j][0] )
			{
				$i += 1;
				next;
			}

			# i.L > j.R
			if ( $aa->{$key}[$i][0] > $bb->{$key}[$j][1] )
			{
				$j += 1;
				next;
			}

			# i.L <= j.R  and i.R >= j.L
			if (( $aa->{$key}[$i][0] <= $bb->{$key}[$j][1] ) && ( $aa->{$key}[$i][1] >= $bb->{$key}[$j][0] ))
			{
				if ( $in_bb )
				{
					push @arr, $bb->{$key}[$j][2];
				}
				else
				{
					push @arr, [ @{$aa->{$key}[$i]} ];

					$hash{ $key ."_". $aa->{$key}[$i][0] ."_". $aa->{$key}[$i][1] } += 1;
				}

				$j += 1;
				next;
			}

			die "error in logic\n";
		}
	}

	my $size_before = 0;
	my $size_after = 0;
	my $multi_removal = 0;

	if ( $in_bb )
	{
		$size_before = scalar ( keys %{$ref} );

		foreach my $key (@arr)
		{
			delete $ref->{$key};
		}

		$size_after = scalar ( keys %{$ref} );
	}
	else
	{
		foreach my $key ( keys %{$aa} )
		{
			$size_before += scalar ( @{$aa->{$key}} );

			my @new_arr = ();

			foreach my $value ( @{$aa->{$key}} )
			{
				my $test_key = $key ."_". $value->[0] ."_". $value->[1];

				if ( ! exists $hash{$test_key} )
				{
					push @new_arr, [ $value->[0], $value->[1] ];
				}
			}

			@{$aa->{$key}} = @new_arr;

			$size_after += scalar (@new_arr);
		}

		foreach my $key ( keys %hash )
		{
			$multi_removal += ( $hash{$key} - 1 );
		}
	}

	if ($v)
	{
		my $size = scalar @arr;
		print "$label $size\n";
		print "# size before, after and multi removal: $size_before $size_after $multi_removal\n";	
#		print Dumper(\@arr);
	}
}
# ------------------------------------------------
sub HashKeysToArrIntervals
{
	my $ref = shift;

	my %h = ();

	foreach my $key (keys %{$ref})
	{
#		if ( $key =~ /^(.+)_(\d+)_(\d+)_[+-.]$/ or $key =~ /^(.+)_(\d+)_(\d+)_[+-.]_[012.]$/ )
#		{
#			push @{$h{$1}}, [$2, $3, $key];
#		}
#		elsif ( $key =~ /^(.+?)_(\d+).*_(\d+)_[+-.] $/ or $key =~ /^(.+?)_(\d+).*_(\d+)_[+-.]_[012.] $/ )

		if ( $key =~ /^(.+?)_(\d+).*_(\d+)_[+-.]$/ or $key =~ /^(.+?)_(\d+).*_(\d+)_[+-.]_[012.]$/ )
		{
			push @{$h{$1}}, [$2, $3, $key];
		}
		else
			{ die "error, unexpected key format found: \"$key\"\n"; }
	}

	foreach my $key (keys %h)
	{
		my @arr_a = @{$h{$key}};
		my @arr_b = sort{ $a->[0] <=> $b->[0] } (@arr_a);

		$h{$key} = \@arr_b;
	}

	return %h;
}
# ------------------------------------------------
sub Compare3
{
	my $c1 = scalar keys %h1;
	my $c2 = scalar keys %h2;
	my $c3 = scalar keys %h3;
	
	my $match = 0;
	
	my $c1_c2 = 0;
	my $c1_c3 = 0;
	my $c2_c3 = 0;
	
	my $c1_uniq = 0;
	my $c2_uniq = 0;
	my $c3_uniq = 0;
	

	foreach my $key ( keys %h1 )
	{
		if ( exists $h2{$key} and exists $h3{$key} )
		{
			$match += 1;
			$shared_123{$key} += 1 if $out_file;
		}
		elsif ( exists $h2{$key} and  ! exists $h3{$key} )
		{
			$c1_c2 += 1;
			$shared_12{$key} += 1 if $out_file;
		}
		elsif (! exists $h2{$key} and  exists $h3{$key} )
		{
			$c1_c3 += 1;
			$shared_13{$key} += 1 if $out_file;
		}
		else
		{
			$unique_1{$key} += 1 if $out_file;
		}
	}

	foreach my $key ( keys %h2 )
	{
		if ( ! exists $h1{$key} and exists $h3{$key} )
		{
			$c2_c3 += 1;
			$shared_23{$key} += 1 if $out_file;
		}
		
		if ($out_file)
		{
			if ( ! exists $h1{$key} and ! exists $h3{$key} )
			{
				$unique_2{$key} += 1 if $out_file;
			}
		}
	}

	if ( $out_file )
	{
		foreach my $key ( keys %h3 )
                {
			if ( ! exists $h1{$key} and ! exists $h2{$key} )
			{
				$unique_3{ $key } += 1;
			}
		}
	}

	$c1_uniq = $c1 - $match - $c1_c2 - $c1_c3;
        $c2_uniq = $c2 - $match - $c1_c2 - $c2_c3;
        $c3_uniq = $c3 - $match - $c1_c3 - $c2_c3;

	die "error, no coding intervals in file: $f1\n" if (!$c1);
	die "error, no coding intervals in file: $f2\n" if (!$c2);
	die "error, no coding intervals in file: $f3\n" if (!$c3);

	if ($v)
	{
		print "\n";
		print "#in\tmatch\tunique\t\%match\tCDS\n"  if $compare_cds;
		print "#in\tmatch\tunique\t\%\tIntron\n"    if $compare_introns;
		print "#in\tmatch\tunique\t\%\tStarts\n"    if $compare_starts;
		print "#in\tmatch\tunique\t\%\tStops\n"     if $compare_stops;
		print "#in\tmatch\tunique\t\%\tDonors\n"    if $compare_donors;
		print "#in\tmatch\tunique\t\%\tAcceptors\n" if $compare_acceptors;
		print "#in\tmatch\tunique\t\%\tSingle-CDS\n" if $compare_single;
		print "#in\tmatch\tunique\t\%\tInitial-CDS\n" if $compare_initial;
		print "#in\tmatch\tunique\t\%\tInternal-CDS\n" if $compare_internal;
		print "#in\tmatch\tunique\t\%\tTerminal-CDS\n" if $compare_terminal;
		print "#in\tmatch\tunique\t\%\tTranscript-CDS\n" if $compare_trans;
		print "#in\tmatch\tunique\t\%\tGene-CDS\n" if $compare_gene;
	}

	print "\n";
	print $c1 ."\t". $match ."\t". "0"    ."\t". $c1_c2 ."\t". $c1_c3 ."\t". $c1_uniq ."\t". sprintf( "%.2f", 100.0*$match/$c1 ) ."\t". $f1 ."\n";
	print $c2 ."\t". $match ."\t". $c1_c2 ."\t". "0"    ."\t". $c2_c3 ."\t". $c2_uniq ."\t". sprintf( "%.2f", 100.0*$match/$c2 ) ."\t". $f2 ."\n";
	print $c3 ."\t". $match ."\t". $c1_c3 ."\t". $c2_c3 ."\t". "0"    ."\t". $c3_uniq ."\t". sprintf( "%.2f", 100.0*$match/$c3 ) ."\t". $f3 ."\n";
	print "\n";

	if ($out_file)
	{
		PrintKeys( \%unique_1,  "unique 1",   $out_file ) if $unique1;
		PrintKeys( \%unique_2,  "unique 2",   $out_file ) if $unique2;
		PrintKeys( \%unique_3,  "unique 3",   $out_file ) if $unique3;
		PrintKeys( \%shared_12, "shared 1-2", $out_file ) if $shared12;
		PrintKeys( \%shared_13, "shared 1-3", $out_file ) if $shared13;
		PrintKeys( \%shared_23, "shared 2-3", $out_file ) if $shared23;
		PrintKeys( \%shared_123, "shared 1-2-3", $out_file ) if $shared123;

		if ($debug)
		{
			TestCount( $match, \%shared_123 );
			TestCount( $c1_c2, \%shared_12 );
			TestCount( $c1_c3, \%shared_13 );
			TestCount( $c2_c3, \%shared_23 );
			TestCount( $c1_uniq, \%unique_1 );
			TestCount( $c2_uniq, \%unique_2 );
			TestCount( $c3_uniq, \%unique_3 );
			TestCount( $c1, \%shared_123, \%shared_12, \%shared_13, \%unique_1 );
			TestCount( $c2, \%shared_123, \%shared_12, \%shared_23, \%unique_2 );
			TestCount( $c3, \%shared_123, \%shared_13, \%shared_23, \%unique_3 );
        	}
	}
}
# ------------------------------------------------
sub PrintKeys
{
	my $ref = shift;
	my $label = shift;
	my $name = shift;

	open( my $OUT, ">", $name ) or die;

	print $OUT "# $label\n";

	foreach my $key (sort keys %{$ref})
	{
		if ( ! $original )
		{
			print $OUT "$key\n";
		}
		else
		{
			my $txt = PrintOriginal( $key, $original );
			if ( $txt !~ /\n$/ )
			{
				$txt .= "\n";
			}
			print $OUT $txt;
		}
	}

	close $OUT;
}
# ------------------------------------------------
sub PrintOriginal
{
	my $key = shift;
	my $id = shift;

	if ( $id == 1 )
	{
		return $h1{$key};
	}
	elsif ( $id == 2 )
	{
		return $h2{$key};
	}
	elsif ( $id == 3 )
	{
		return $h3{$key};
	}
	else
		{ die "error, unexpected value found in PrintOriginal\n"; }
}
# ------------------------------------------------
sub Compare2
{
	my $total_c1 = scalar keys %h1;
	my $total_c2 = scalar keys %h2;
	my $match = 0;
	my $c1_uniq = 0;
	my $c2_uniq = 0;
	
	foreach my $key ( keys %h1 )
	{
		if ( exists $h2{$key} )
		{
			$match += 1;
			$shared_12{$key} += 1 if $out_file;
		}
		else
		{
			$unique_1{$key} += 1 if $out_file
		}
	}

	if ( $out_file )
	{
		foreach my $key ( keys %h2 )
		{
			if ( ! exists $h1{$key} )
			{
				$unique_2{$key} += 1;
			}
		}
	}

	$c1_uniq = $total_c1 - $match;
	$c2_uniq = $total_c2 - $match;

	die "error, no intervals in file: $f1\n" if (!$total_c1);
	die "error, no intervals in file: $f2\n" if (!$total_c2);

	if ($v)
	{
		print "\n";
		print "#in\tmatch\tunique\t\%match\tCDS\n"  if $compare_cds;
		print "#in\tmatch\tunique\t\%\tIntron\n"    if $compare_introns;
		print "#in\tmatch\tunique\t\%\tStarts\n"    if $compare_starts;
		print "#in\tmatch\tunique\t\%\tStops\n"     if $compare_stops;
		print "#in\tmatch\tunique\t\%\tDonors\n"    if $compare_donors;
		print "#in\tmatch\tunique\t\%\tAcceptors\n" if $compare_acceptors;
		print "#in\tmatch\tunique\t\%\tSingle-CDS\n" if $compare_single;
		print "#in\tmatch\tunique\t\%\tInitial-CDS\n" if $compare_initial;
		print "#in\tmatch\tunique\t\%\tInternal-CDS\n" if $compare_internal;
		print "#in\tmatch\tunique\t\%\tTerminal-CDS\n" if $compare_terminal;
		print "#in\tmatch\tunique\t\%\tTranscript-CDS\n" if $compare_trans;
		print "#in\tmatch\tunique\t\%\tGene-CDS\n" if $compare_gene;
	}

	print "\n";
	print  $total_c1 ."\t". $match ."\t". $c1_uniq ."\t". sprintf( "%.2f", 100.0*$match/$total_c1 ) ."\t". $f1 ."\n";
	print  $total_c2 ."\t". $match ."\t". $c2_uniq ."\t". sprintf( "%.2f", 100.0*$match/$total_c2 ) ."\t". $f2 ."\n";
	print "\n";

	if ( $out_file )
	{
		PrintKeys( \%unique_1,  "unique 1",   $out_file ) if $unique1;
		PrintKeys( \%unique_2,  "unique 2",   $out_file ) if $unique2;
		PrintKeys( \%shared_12, "shared 1-2", $out_file ) if $shared12;

		if ($debug)
		{
			TestCount( $match, \%shared_12 );
			TestCount( $c1_uniq, \%unique_1 );
			TestCount( $c2_uniq, \%unique_2 );
			TestCount( $total_c1, \%shared_12, \%unique_1 );
			TestCount( $total_c2, \%shared_12, \%unique_2 );
		}
	}
}
# ------------------------------------------------
sub TestCount
{
	my $count = shift;
	my $ref1 = shift;
	my $ref2 = shift;
	my $ref3 = shift;
	my $ref4 = shift;

	my $size = scalar (keys %{$ref1});

	$size += scalar (keys %{$ref2}) if ($ref2);
	$size += scalar (keys %{$ref3}) if ($ref3);
	$size += scalar (keys %{$ref4}) if ($ref4);

	die "error in counting: $count $size\n" if ( $size != $count );
}
# ------------------------------------------------
sub LoadIntervals
{
	my ($name, $min) = @_;

	my %h = ();

	my $in_regions = 0;

	open( my $IN, $name ) or die "error on open file $name: $!\n";
	while( my $line = <$IN> )
	{
		next if ( $line =~ /^\s*#/ );
		next if ( $line =~ /^\s*$/ );
		next if ( $line !~ /\t/ );
		
		if( $line =~ /^(\S+)\t\S+\t(\S+)\t(\d+)\t(\d+)\t\S+\t([-+.])\t(\S+)\s*/ )
		{
			my $id     = $1;
			my $type   = $2;
			my $start  = $3;
			my $end    = $4;
			my $strand = $5;
			my $ph     = $6;

			next if ( $end - $start + 1 < $min ); 

			push @{$h{$id}}, [$start, $end];
			$in_regions += 1;
		}
		else
			{ die "error, unexpected line format found: $line"; }
	}
	close $IN;

	my $out_regions = 0;

	foreach my $key (keys %h)
	{
		my @arr_a = @{$h{$key}};
		my @arr_b = sort{ $a->[0] <=> $b->[0] } (@arr_a);
		my @arr_c = ();

		my $size = scalar @arr_b;

		for( my $i = 0; $i < $size; )
		{
			my $L = $arr_b[$i][0];
			my $R = $arr_b[$i][1];

			my $j = 0;

			for( $j = $i+1; $j < $size; $j += 1 )
			{
				if ( $R >= $arr_b[$j][0] )
				{
					if ( $R < $arr_b[$j][1] )
					{
						$R = $arr_b[$j][1];
						print "# merged $i $j $L $R $arr_b[$j][0] $arr_b[$j][1]\n" if $debug
					}
					else
					{
						print "# full overlap found $i $j $L $R $arr_b[$j][0] $arr_b[$j][1]\n" if $debug
					} 
				}
				else
				{
					last;
				}
			}

			$i = $j;

			push @arr_c, [$L, $R];
		}

		$h{$key} = \@arr_c;

		$out_regions += (scalar @arr_c);
	}

	if ($v)
	{
		print "# regions in: $in_regions\n";
		print "# regions after merging: $out_regions\n";

		if ($debug)
		{
			my $sum = 0;
			foreach my $key (keys %h )
			{
				$sum += (scalar @{$h{$key}});
				print "# per sequence: ". $key ."\t". (scalar @{$h{$key}}) ."\n";

				foreach my $value (@{$h{$key}})
				{
					print "###". $key ."\t". $value->[0] ."\t". $value->[1] ."\n";
				}
			}
			print "# total $sum\n";
		}
	}

#	print Dumper(\%h);

	return %h;
}
# ------------------------------------------------
sub ParseGFF
{
	my ($fname, $ref, $tr2g) = @_;

	# $fname - file name; file in GFF-like format
	# $ref - reference on hash; key - created from record; 
	# $tr2g - reference on hash; key - transcript key; value - gene ID

	open( my $IN, $fname ) or die "error on open file $fname: $!\n";
	while( my $line = <$IN> )
	{
		next if ( $line =~ /^\s*#/ );
		next if ( $line =~ /^\s*$/ );
		next if ( $line !~ /\t/ );
		
		if ( $line =~ /^([^\t]+)\t[^\t]+\t([^\t]+)\t(\d+)\t(\d+)\t\S+\t([-+.])\t([0123.])(\t(.*)|\s*)$/ )
		{
			my $id     = $1;
			my $type   = $2;
			my $start  = $3;
			my $end    = $4;
			my $strand = $5;
			my $ph     = $6;
			my $attr   = $7;

			# attr is not required in GFF
			# if not found - set it to empty

			$attr =~ s/^\t// if ( defined $attr );
			$attr = ''       if ( !defined $attr );

			if ( $compare_cds )
			{
				next if ( $type ne "CDS" );
			}
			elsif ( $compare_exon )
			{
				next if ( $type ne "exon" );
			}
			elsif ( $compare_introns )
			{
				if ( $type =~ /^[Ii]ntron$/ ) {;} else {next;}
			}
			elsif ( $compare_starts )
			{
				if ( $type =~ /^[Ss]tart_codon$/ ) {;} else {next;}
			}
			elsif ( $compare_stops )
			{
				if ( $type =~ /^[Ss]top_codon$/ ) {;} else {next;}
			}
			elsif ( $compare_donors )
			{
				if ( $type =~ /^[Ii]ntron$/ ) {;} else {next;}
			}
			elsif ( $compare_acceptors )
			{
				if ( $type =~ /^[Ii]ntron$/ ) {;} else {next;}
			}
			elsif ( $compare_single )
			{
				die "error, CDS whithout label was found: $attr\n" if (( $type eq "CDS" )and not( $attr =~ /cds_type=/ ));
				if (( $type eq "CDS" ) and ( $attr =~ /(cds_type=[Ss]ingle|cds_type \"[Ss]ingle\")/ )) {;} else {next;}
			}
			elsif ( $compare_initial )
			{
				die "error, CDS whithout label was found: $attr\n" if (( $type eq "CDS" )and not( $attr =~ /cds_type=/ ));
				if (( $type eq "CDS" ) and ( $attr =~ /(cds_type=[Ii]nitial|cds_type \"[Ii]nitial\")/ )) {;} else {next;}
			}
			elsif ( $compare_internal )
			{
				die "error, CDS whithout label was found: $attr\n" if (( $type eq "CDS" )and not( $attr =~ /cds_type=/ ));
				if (( $type eq "CDS" ) and ( $attr =~ /(cds_type=[Ii]nternal|cds_type \"[Ii]nternal\")/ )) {;} else {next;}
			}
			elsif ( $compare_terminal )
			{
				die "error, CDS whithout label was found: $attr\n" if (( $type eq "CDS" )and not( $attr =~ /cds_type=/ ));
				if (( $type eq "CDS" ) and ( $attr =~ /(cds_type=[Tt]erminal|cds_type \"[Tt]erminal\")/ )) {;} else {next;}
			}
			elsif ( $compare_trans )
			{
				next if ( $type ne "CDS" );
			}
			elsif ( $compare_gene )
			{
				next if ( $type ne "CDS" );
			}
			elsif ( $compare_multiGene )
			{
				die "error, CDS whithout label was found: $attr\n" if (( $type eq "CDS" )and not( $attr =~ /cds_type=/ ));
				if (( $type eq "CDS") and not ( $attr =~ /(cds_type=[Ss]ingle|cds_type \"[Ss]ingle\")/ )) {;} else {next;}
			}
			elsif ( $compare_singleGene )
			{
				die "error, CDS whithout label was found: $attr\n" if (( $type eq "CDS" )and not( $attr =~ /cds_type=/ ));
				if (( $type eq "CDS") and ( $attr =~ /(cds_type=[Ss]ingle|cds_type \"[Ss]ingle\")/ )) {;} else {next;}
			}
			else
				{ next; }

			if ( $v )
			{
				print "warning, strand is not defined: $line" if ($strand eq ".");
			}

			# start of unique key for most of the comparisons
			# 
			my $key = $id ."_". $start ."_". $end ."_". $strand;
			
			# donors and acceptors have only one coordinate, so restart the key in the new format

			if ( $compare_donors )
			{
				if ( $strand eq "+" )
				{
					$key = $id ."_". $start ."_". $strand;
				}
				elsif ( $strand eq "-" )
				{
					$key = $id ."_". $end   ."_". $strand;
				}
				else {die "error, strand is missing for compare donors from introns\n";}
			}

			if ( $compare_acceptors )
			{
				if ( $strand eq "+" )
				{
					$key = $id ."_". $end   ."_". $strand;
				}
				elsif ( $strand eq "-" )
				{
					$key = $id ."_". $start ."_". $strand;
				}
				else {die "error, strand is missing for compare acceptors from introns\n";}
			}

			if ( ! $no_phase )
			{
				$key .= "_". $ph;
			}

			# insert key->value into hash

			if ( !( $compare_trans or $compare_gene or $compare_multiGene or $compare_singleGene ))
			{
				if ( ! $original )
				{
					$ref->{$key} += 1;
				}
				else
				{
					$ref->{$key} .= $line;
				}
			}

			# Gene and transcript IDs are parsed only from the GTF format
			# This procedure creates a temporary key-to-true-key hash
			# Hash should be reversed after processing of records from the input file
			# This is done at the end of the sub

			if ( $compare_trans or $compare_gene or $compare_multiGene or $compare_singleGene )
			{
				my $gene_id = '';
				my $trans_id = '';

				if ( $line =~ /\sgene_id \"(\S+)\"\;/ )
				{
					$gene_id = $1;
				}

				if ( $line =~ /\stranscript_id \"(\S+)\"\;/ )
				{
					$trans_id = $1;
				}

				if ( !$gene_id or !$trans_id )
					{ die "error, GTF formatted file with gene and transcript ID is expected: $line\n"; }

				$ref->{$trans_id} .= ($key ."\t");
				$tr2g->{$trans_id} = $gene_id;
			}
		}
		else
		{
			print "warning, unxpected line format found in gff: $line\n";
		}
	}	
	close $IN;

	if ($v)
	{
		print "# CDS in file $fname: ".          (scalar keys %$ref) ."\n" if $compare_cds;
		print "# Introns in file $fname: ".      (scalar keys %$ref) ."\n" if $compare_introns;
		print "# Starts in file $fname: ".       (scalar keys %$ref) ."\n" if $compare_starts;
		print "# Stops in file $fname: ".        (scalar keys %$ref) ."\n" if $compare_stops;
		print "# Donors in file $fname: ".       (scalar keys %$ref) ."\n" if $compare_donors;
		print "# Acceptors in file $fname: ".    (scalar keys %$ref) ."\n" if $compare_acceptors;
		print "# Single-CDS in file $fname: ".   (scalar keys %$ref) ."\n" if $compare_single;
		print "# Initial-CDS in file $fname: ".  (scalar keys %$ref) ."\n" if $compare_initial;
		print "# Internal-CDS in file $fname: ". (scalar keys %$ref) ."\n" if $compare_internal;
		print "# Terminal-CDS in file $fname: ". (scalar keys %$ref) ."\n" if $compare_terminal;
		print "# Exons in file $fname: ".        (scalar keys %$ref) ."\n" if $compare_exon;
	}

	if ( $compare_trans or $compare_gene or $compare_multiGene or $compare_singleGene )
	{
		if ($v)
		{
			print "# Transcripts-CDS in file $fname: ". (scalar keys %$ref) ."\n" if ( $compare_trans or $compare_gene );
			print "# Transcripts-multi-CDS in file $fname: ". (scalar keys %$ref) ."\n" if $compare_multiGene;
			print "# Transcripts-single-CDS in file $fname: ". (scalar keys %$ref) ."\n" if $compare_singleGene;
		}

		SortValueWithTabSep($ref);
		%{$ref} = ReverseKeyValueWithRevInfo( $ref );
		SortValueWithTabSep($ref);

		CheckForSharedTranscripts($ref, $tr2g);

		if ($v)
		{
			print "## After duplication removal\n";
			print "## Transcripts-CDS in file $fname: ". (scalar keys %$ref) ."\n" if ( $compare_trans or $compare_gene );
			print "## Transcripts-multi-CDS in file $fname: ". (scalar keys %$ref) ."\n" if $compare_multiGene;
			print "## Transcripts-single-CDS in file $fname: ". (scalar keys %$ref) ."\n" if $compare_singleGene;
		}
	}
}
# ------------------------------------------------
sub CheckForSharedTranscripts
{
	my $ref = shift;
	my $tr2g = shift;

	my $shared = 0;

	foreach my $key (keys %{$ref})
	{
		my %gid = ();
		my $size = 0;

		# tid's are in arr
		my @arr = split( "\t", $ref->{$key} );

		foreach my $id (@arr)
		{
			$gid{ $tr2g->{$id} } .= $id ."\t";
		}

		$size = scalar (keys %gid);

		if ( $size != 1 )
		{
			print "# warning, transcript $key found in $size genes\n" if $debug;
			print Dumper(\%gid) if $debug;
			$shared += $size;
		}
	}

	print "## Shared transcripts: $shared\n" if $v;
}
# ------------------------------------------------
sub SortValueWithTabSep
{
	my $ref = shift;

	foreach my $key (keys %{$ref})
	{
		my @arr = split( "\t", $ref->{$key} );
		my %tmp = ();
		foreach my $value (@arr)
		{
			$tmp{$value} += 1;
		}
		@arr = sort keys %tmp;
		$ref->{$key} = join( "\t", @arr );
	}
}
# ------------------------------------------------
sub ReverseKeyValueWithRevInfo
{
	my $ref = shift;

	my %h = ();

	foreach my $key (sort keys %{$ref})
	{
		if ( !exists $h{ $ref->{$key} } )
		{
			$h{ $ref->{$key} } = $key;
		}
		else
		{
			$h{ $ref->{$key} } .= ("\t". $key);
		}
	}

	return %h;
}
# ------------------------------------------------
sub CheckBeforeRun
{
	die "error, file not found: option --f1 $f1\n" if( ! -e $f1 );
	die "error, file not found: option --f2 $f2\n" if( ! -e $f2 );

	die "error, file is missing: option --f3\n" if( $f3 and ( ! -e $f3 ));

	if ( $out_file )
	{
		die "error, output file name matches input file\n" if (( $out_file eq $f1 ) or ( $out_file eq $f1 ));
		die "error, output file name matches input file\n" if ( $f3 and ( $out_file eq $f3 ));
	}

	die "error, file not found: option --pseudo $pseudo\n" if ( $pseudo and ( ! -e $pseudo ));
	die "error, file not found: option --mask $masked\n" if ( $masked and ( ! -e $masked ));

	if ( $original )
	{
		if ( $original == 1 or $original == 2 )
			{;}
		elsif ( $f3 and $original == 3 )
			{;}
		else
			{ die "error, value of --original parameter is out of allowed range: $original\n"; }
	}
}
# ------------------------------------------------
sub ParseCMD
{
	my $opt_results = GetOptions
	(
		'f1=s'   => \$f1,
		'f2=s'   => \$f2,
		'f3=s'   => \$f3,
		'out=s'  => \$out_file,
		'cds'       => \$compare_cds,
		'introns'   => \$compare_introns,
		'don'       => \$compare_donors,
		'acc'       => \$compare_acceptors,
		'starts'    => \$compare_starts,
		'stops'     => \$compare_stops,
		'no_phase'  => \$no_phase,
		'verbose'   => \$v,
		'debug'     => \$debug,
		'shared123' => \$shared123,
		'shared12'  => \$shared12,
		'shared13'  => \$shared13,
		'shared23'  => \$shared23,
		'unique1'   => \$unique1,
		'unique2'   => \$unique2,
		'unique3'   => \$unique3,
		'original=i' => \$original,
		'pseudo=s'   => \$pseudo,
		'masked=s'   => \$masked,
		'min_mask=i' => \$min_masked,
		'single'     => \$compare_single,
		'initial'    => \$compare_initial,
		'internal'   => \$compare_internal,
		'terminal'   => \$compare_terminal,
		'gene'       => \$compare_gene,
		'multiGene'  => \$compare_multiGene,
		'singleGene' => \$compare_singleGene,
		'trans'      => \$compare_trans,
		'exon'       => \$compare_exon,
	);

	die "error on command line\n" if( !$opt_results );
	die "error, unexpected argument found on command line\n" if( @ARGV > 0 );

	my $count = 0;
	$count += 1 if $compare_cds;
	$count += 1 if $compare_introns;
	$count += 1 if $compare_donors;
	$count += 1 if $compare_acceptors;
	$count += 1 if $compare_starts;
	$count += 1 if $compare_stops;
	$count += 1 if $compare_single;
	$count += 1 if $compare_initial;
	$count += 1 if $compare_internal;
	$count += 1 if $compare_terminal;
	$count += 1 if $compare_gene;
	$count += 1 if $compare_multiGene;
	$count += 1 if $compare_singleGene;
	$count += 1 if $compare_trans;
	$count += 1 if $compare_exon;

	die "error, more than one comparison type was specified on the command line: $count\n" if ($count > 1);

	# set default comparison mode
	$compare_cds = 1 if ($count == 0);

	$v = 1 if $debug;
}
# ------------------------------------------------
sub Usage
{
	print qq(
Usage:
$0  --f1 [name]  --f2 [name]

Compare intervals from input files:

   --f1  [name]  file in GFF-like format
   --f2  [name]  file in GFF-like format

   By default, a comparison is done for records with CDS features (coding part of the exon).

Optional:
   --f3  [name]  file in GFF-like format

   --pseudo [name] file with pseudogene coordinates in GFF-like format;
                   with this option program:
                     a. removes 'pseudo' label from pseudo regions overlapping with --f1
                     b. excludes from comparison regions from --f2 and --f3 that overlap pseudo regions

   --masked [name] file with coordinates of masked regions;
                   with this option any region that overlaps masking
                   is excluded from the comparison
   --min_mask [$min_masked] minimum length of the masked region to use

Default comparison is done for the 'CDS' type

   --cds         compare fields of 'CDS' type
   --introns     compare fields of 'intron' type
   --starts      compare fields of 'start_codon' type
   --stops       compare fields of 'stop_codon' type
   --don         compare donors using the 'intron' type
   --acc         compare acceptors using the 'intron' type
   --exon        compare fields of 'exon' type

   Using optional 'cds_type' label in attribute field of GFF-like format

   --single      compare fields of 'CDS' type with 'single' label
   --initial     compare fields of 'CDS' type with 'initial' label
   --internal    compare fields of 'CDS' type with 'internal' label
   --terminal    compare fields of 'CDS' type with 'terminal' label

   Using gene and transcript IDs from the GTF formatted file

   --trans       compare transcripts (CDS bases)
   --gene        compare genes (CDS based)
                 gene is called 'found' if at least one transcript was matched exactly

   --multigene   compare only multi-CDS genes; label based
   --singlegene  compare only single-CDS genes; label based

   --no_phase    ignore phase of record in comparison

Save GFF-like record keys into the output file based on the record status:

   --out [name]  output file name

   --shared12    output shared by 1 and 2, not 3
   --shared13    output shared by 1 and 3, not 2
   --shared23    output shared by 2 and 3, not 1

   --shared123   output shared by 1-2-3

   --unique1     output unique in 1
   --unique2     output unique in 2
   --unique3     output unique in 3

   --original [number]  1 or 2 of 3; output full record (instead of key) as in original input files
                 1 corresponds to the file with --f1 name, etc.
General:
   --verbose
   --debug

Developer:
   --

Version: $VERSION

);
	exit 1;
}
# ------------------------------------------------

