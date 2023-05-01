#!/usr/bin/perl

use strict;
use warnings;

use lib '../';

use OriginCoord;

print STDERR "Loading Origin coordinates\n";
my $origins = OriginCoord->new();

my %chromosomes = $origins->get_chromosomes();
my %oristart = $origins->get_origin_starts();
my %oriend = $origins->get_origin_ends();


# get chromosome length
my %chrlen = ();
my $chrfile = "../chrom_len_saccer3.txt";
open (CHROM, $chrfile ) || die "couldn't open file\n";
while ( my $ch = <CHROM> )
{
        chomp $ch;
        my @val = split /\t/, $ch;
        $chrlen{$val[0]} = $val[1];
}

my $chr = "";
my %strandlookup;
my $mitomutcount = 0;

my $bedfile = "saccer3_origins_8000bpflank.bed";
open ( BED, ">$bedfile") || die "couldn't open file $bedfile\n";

# distance from origin to count
my $window = 8000;

foreach my $chr ( sort keys %chromosomes )
{
	my $oricount = scalar @{$chromosomes{$chr}};
        print STDERR "Starting to process $chr, which has $oricount origins\n";
	%strandlookup = ();
	my $skipped = 0;
	foreach my $acc ( @{$chromosomes{$chr}} )
	{
                if ( not exists $oristart{$acc} || not exists $oriend {$acc} )
                {
			$skipped++;
                        next;
               	}

		# indicate strand of leading strand --> left of origin, + strand is leading
		my $leftstrand = "+";
		my $leftstart = $oristart{$acc} - $window;
		$leftstart--; # bed formatting
		if ( $leftstart < 0 )
		{
			$leftstart = 0;
		}

		my $leftend = $oristart{$acc} - 1;
		# strand indicates leading strand
		print BED "$chr\t$leftstart\t$leftend\t$acc\t.\t$leftstrand\n";

		# indicate strand of leading strand --> right of origin, - strand is leading
		my $rightstrand = "-";
		my $rightstart = $oriend{$acc}; # because of bed formatting

		my $rightend = $oriend{$acc} + $window;
		if ( $rightend > $chrlen{$chr} )
		{
			$rightend = $chrlen{$chr};
		}
		print BED "$chr\t$rightstart\t$rightend\t$acc\t.\t$rightstrand\n";

	}

	print STDERR "$skipped skipped origins\n";
}

system ( "printf 'saccer3_origins_8000bpflank.bed\n' | perl ../split_strands.pl" ); 

system ("perl ../merge_genecoords.pl <saccer3_origins_8000bpflank_plusstrand.bed >saccer3_origins_8000bpflank_plusmerge.bed");
system ("perl ../merge_genecoords.pl <saccer3_origins_8000bpflank_minusstrand.bed >saccer3_origins_8000bpflank_minusmerge.bed");

system ( "cat saccer3_origins_8000bpflank_plusmerge.bed saccer3_origins_8000bpflank_minusmerge.bed >saccer3_origins_8000bpflank_bothmerge.bed");
system ( "cat saccer3_origins_8000bpflank_bothmerge.bed | sort -k1,1 -k2,2n -k3,3n >saccer3_origins_8000bpflank_bothmerge_sort.bed" );

my $mergefile = $bedfile;
$mergefile =~ s/\.bed/_noverlap\.bed/;

# subtract overlapping regions

system ( "perl ../removerlap_genecoord.pl <saccer3_origins_8000bpflank_bothmerge_sort.bed >$mergefile" );

my $fa_file = $mergefile;
$fa_file =~ s/\.bed/\.fa/;
system ("fastaFromBed -s -name -fi ../saccer3_genome.fa -bed $mergefile -fo $fa_file" );

# process sequence in fasta file
open ( FASTA, $fa_file ) || die "Couldn't open fasta file\n";

my %trinucount;
my $genelen = 0;
my %lagging;
my %leading;
while ( my $line = <FASTA> )
{
        chomp $line;
	if ( $line =~ /::chr[XIVM]+:([0-9]+)-([0-9]+)/ )
	{
		$genelen = abs( $1 - $2 );
		next;	# header of fasta
	}

	if ( length($line) != ( $genelen ) )
	{
		die "FASTA line of wrong length (should be $genelen): $line\n";
	}
	
	for ( my $i = 0; $i < (length($line) - 2); $i++ )
	{
		my $triseq = substr( $line, $i, 3 );
		if ( $triseq =~ /^[ATGC][ATGC][ATGC]$/ )
		{
			$trinucount{$triseq}++;
		}
		else
		{
			die "Error with $triseq\n";
		}

		my $mid = substr $triseq, 1, 1;
	
		if ( $mid =~ /[CT]/ )
		{
			$leading{$triseq}++;
		}
		elsif ( $mid =~ /[AG]/ )
		{
			$triseq = reverse $triseq;
			$triseq =~ tr/ACGT/TGCA/;
			$lagging{$triseq}++;
		}
		else
		{
			die "Error if double mutant: $triseq\n";
		}
	}
}

my $outfile = $mergefile;
$outfile =~ s/\.bed/_trinucount\.txt/;
open (OUT, ">$outfile" ) || die "Couldn't open OUT file: $outfile\n";

print OUT "Trinucleotide sequence\tCount in Iyer genes\n";
foreach my $tri ( sort keys %trinucount )
{
	print OUT "$tri\t$trinucount{$tri}\n";
}
close ( OUT );

my $sigfile = $bedfile;
$sigfile =~ s/\.bed/_strand_trinucfreq\.txt/;
open (SIG, ">$sigfile" ) || die "Couldn't open OUT file: $sigfile\n";
print SIG "Trinucleotide sequence\tLagging count\tLeading count\n";
foreach my $tri ( sort keys %lagging )
{
	print SIG "$tri\t$lagging{$tri}\t$leading{$tri}\n";
}
