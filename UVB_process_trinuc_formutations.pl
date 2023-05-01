#!/usr/bin/perl -w

use strict;
use warnings;

my $header = <STDIN>;
print STDERR "Header: $header\n";

my $trifile = "trinucontext_UVB_mutations_yeast.bed";
open ( TRI, ">$trifile" ) || die "couldn't open $trifile\n";
my $fafile = "trinucontext_UVB_mutations_yeast.fa";

while( my $line = <STDIN> )
{
	chomp $line;
	if ( $line eq "")
	{
		next;
	}
	my @fields = split /\t/, $line;

	if ( $fields[11] =~ /SNV/ )
	{
        	my $chr = $fields[2];
		if ( $chr eq "chrMito" )
		{
			$chr = "chrM";
		}
	        my $mut_end = $fields[4];

        	# convert to 0-based bedfile format from 1-based mutation call
	        my $mut_start = $mut_end - 1;

		if ( $fields[3] != $mut_end )
		{
			die "Error in mutstart and mutend!\n";
		}
		# Get strand and base information for mutation
		my $origbase = $fields[12];
		my $id = $fields[8];
		my $newbase = $fields[13]; 
		my $strand = "";
		my $refbase = "";
		my $mutbase = "";
		if ( $origbase eq $mutbase )
		{
			die "Reference base equals mutbase: $line\n";
		}
		if ( $origbase eq "T" || $origbase eq "C" )
		{
			$strand = "+";
			$refbase = $origbase;
			$mutbase = $newbase;
		}
		elsif ( $origbase eq "A" || $origbase eq "G" )
		{
			$strand = "-";
			$refbase = $origbase;
			$refbase =~ tr/ATGC/TACG/;
			$mutbase = $newbase;
			$mutbase =~ tr/ATGC/TACG/;
		}
		else
		{
			print STDERR "ERROR: reference base is $origbase\n";
			next;
		}

		# expand to trinuc context
		$mut_start -= 1;
		$mut_end += 1;
        	print TRI "${chr}\t$mut_start\t$mut_end\t$id\t$mutbase\t$strand\n";
	}
}
system ("fastaFromBed -s -fi ../saccer3_genome.fa -bed $trifile -fo $fafile");
