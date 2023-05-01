#!/usr/bin/perl

use strict;
use warnings;

my %trimuts;
my %trimutsub;
my $totalmutcount = 0;

print STDERR "Please enter the name of mutation bed file\n";
my $mutfile = <STDIN>;
chomp $mutfile;

my $outfile = $mutfile;
$outfile =~ s/\.bed/_trinucmutation_strand_background\.txt/ || die "Mutation file isn't bed format\n";
open ( OUT, ">$outfile" ) || die "Couldn't open OUT file\n";

my $sigfile = $mutfile;
$sigfile =~ s/\.bed/_trinucsignatures.txt/;
open ( SIG, ">$sigfile" );

open ( MUT, $mutfile ) || die "Couldn't open mutation bed file\n";

print STDERR "Please enter the name of the trinucleotide count file for genome\n";
my $trifile = <STDIN>;
chomp $trifile;

open ( TRINUC, $trifile ) || die "Couldn't open trinuc background file\n";

my %trinucs;
# get rid of headers;
my $head1 = <TRINUC>;
my $head2 = <TRINUC>;
my $head3 = <TRINUC>;
my $tricount = 0;
my $totalnucs = 0;
my $mitocount = 0;
while ( my $line = <TRINUC> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	if ( $fields[0] =~ /^([ATGCN]{3})$/ )
	{
		my $context = $1;
		$trinucs{$context} = $fields[1];
		$totalnucs += $fields[1];
		$tricount++;
	}
	else
	{
		die "Misformatted line: $line\n";
	}
}

while ( my $line = <MUT> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	if ( $fields[3] =~ /^([ATGC]{3})$/ )
	{
		my $tri = $1;
		if ( $fields[0] =~ /chr[IXV]+/ )
		{
			$trimuts{$tri}++;
			$trimutsub{$tri}{$fields[4]}++;
			$totalmutcount++;
		}
		elsif ( $fields[0] =~ /chrM/ )
		{
			$mitocount++;
		}
		else
		{
			die "Weird chromosome: $fields[0]\n";
		}
	}
	else
	{
		die "following line doesn't have trinuc format: $line\n";

	}
}

print STDERR "Mitochrondrial mutations: $mitocount\n";

my $overaldensity = 1.0 * $totalmutcount / $totalnucs;
print OUT "Total count of mutations = $totalmutcount\tTotal count of nucleotides = $totalnucs\toveral mutation density = $overaldensity\tMutation file: $mutfile\n";

print OUT "Trinucleotide context\tMutation count\tTrinuc Count\tFrequency of mutation\n";
foreach my $tri ( sort keys %trimuts )
{
	my $fraction = 1.0 * $trimuts{$tri} / $trinucs{$tri};
	print OUT "$tri\t$trimuts{$tri}\t$trinucs{$tri}\t$fraction\n";
}

print SIG "Total count of mutations = $totalmutcount\tTotal count of nucleotides = $totalnucs\toveral mutation density = $overaldensity\tMutation file: $mutfile\\n";
print SIG "\nTrinucleotide context\tMutant base\tMutation substitution\tMutation count\tTrinuc Count\tFrequency of mutation\n";

foreach my $tri ( sort keys %trimutsub )
{
	foreach my $mutbase ( sort keys %{$trimutsub{$tri}} )
	{
		my $refbase = substr $tri, 1, 1;
		my $sub = $refbase . ">" . $mutbase;
		my $fraction = 1.0 * $trimutsub{$tri}{$mutbase} / $trinucs{$tri};
		print SIG "$tri\t$mutbase\t$sub\t$trimutsub{$tri}{$mutbase}\t$trinucs{$tri}\t$fraction\n";

	}
}
