#!/usr/bin/perl

use strict;
use warnings;

use lib '../';

use OriginCoord;

# ask for probe filename to analyze
print STDERR "Enter filename of sorted mutation bedfile\n";
my $mutbedfile = <STDIN>;
chomp($mutbedfile);

my $outputfile = $mutbedfile;
$outputfile =~ s/\.bed/_ReplicationStrandBins\.txt/ || die "Wrong file type!\n";

open (MUT, "$mutbedfile" ) || die "Couldn't open file: $mutbedfile\n";

print STDERR "Loading Origin coordinates\n";
my $origins = OriginCoord->new();

my %lagbin = ();
my %leadbin = ();
my %intergenicbin = ();

my %trinucref = ();
my %subref = ();

my %chromosomes = $origins->get_chromosomes();
my %oristart = $origins->get_origin_starts();
my %oriend = $origins->get_origin_ends();

my $chr = "";
my %strandlookup;
my $mitomutcount = 0;

# distance from origin to count
my $window = 8000;

while ( my $line = <MUT> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	if ( $fields[0] =~ /(chr[IXV]+)/ )
	{
		my $temp = $1;
		if ( $temp ne $chr )
		{
			$chr = $temp;
			my $oricount = scalar @{$chromosomes{$chr}};
                        print STDERR "Starting to process $temp, which has $oricount origins\n";
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
				for ( my $i = $oristart{$acc} - $window; $i < $oristart{$acc}; $i++ )
				{
					if ( exists $strandlookup{$i} )
					{
						if ( $strandlookup{$i} eq "+" || $strandlookup{$i} eq "-" || $strandlookup{$i} eq "AMBIG" )
						{
							if ( $strandlookup{$i} ne $leftstrand )
							{
								$strandlookup{$i} = "AMBIG";
							}
						}
						else
						{
							die "Weird gene lookup: $strandlookup{$i}\n";
						}
					}
					else
					{	
						$strandlookup{$i} = $leftstrand;
					}
				}

				# indicate strand of leading strand --> right of origin, - strand is leading
				my $rightstrand = "-";
                                for ( my $i = $oriend{$acc} + $window; $i > $oriend{$acc}; $i-- )
                                {
                                        if ( exists $strandlookup{$i} )
                                        {
                                                if ( $strandlookup{$i} eq "+" || $strandlookup{$i} eq "-" || $strandlookup{$i} eq "AMBIG" )
                                                {
                                                        if ( $strandlookup{$i} ne $rightstrand )
                                                        {
                                                                $strandlookup{$i} = "AMBIG";
                                                        }
                                                }
                                                else
                                                {
                                                        die "Weird gene lookup: $strandlookup{$i}\n";
                                                }
                                        }
                                        else
                                        {
                                                $strandlookup{$i} = $rightstrand;
                                        }
                                }   			
			}

			print STDERR "$skipped skipped origins\n";
		}
		my $mutpos = $fields[2];  #use 1-based end coord;
		if ( $mutpos != $fields[1] + 1 )
		{
			die "Error in mut coordinates for line: $line\n";
		}
		my $trinuc = $fields[3];
		
		my $mutbase = $fields[4];
		my $mutstrand = $fields[5];
		if ( $mutstrand ne "+" && $mutstrand ne "-" )
		{
			die "Error in mutation strand in line: $line\n";
		}
		my $refbase = substr $trinuc, 1, 1; 
		my $substitution = $refbase . ">" . $mutbase;

		$trinucref{$trinuc} = 1;
		$subref{$substitution} = 1;
		if ( exists $strandlookup{$mutpos} )
		{			
			if ( $strandlookup{$mutpos} ne "+" && $strandlookup{$mutpos} ne "-" )
			{
				if ( $strandlookup{$mutpos} eq "AMBIG" )	
				{
					$intergenicbin{$trinuc}{$substitution}++;
				}
				else
				{
					die "Error with gene strand\n";
				}
			}
			elsif ( $strandlookup{$mutpos} eq $mutstrand )
			{
				$leadbin{$trinuc}{$substitution}++;
			}	
			else
			{
				$lagbin{$trinuc}{$substitution}++;					
			}
		}
		else
		{
			$intergenicbin{$trinuc}{$substitution}++;
		}
	}
	elsif ( $fields[0] =~ /chrM/ )
	{
		$mitomutcount++;
	}
	else
	{
		die "Error in bed file line: $line\n";
	}
}

open (OUT, ">$outputfile") || die "couldn't open file\n";
#print header
print OUT "Data from file: $mutbedfile\tWindow next to origin = $window\n";

print OUT "Trinucleotide context\tMutant base\tMutation substitution\tLagging strand count\tLeading strand count\tIntergenic count\n";
foreach my $tri ( sort keys %trinucref )
{
	my $trirefbase = substr $tri, 1, 1;
	foreach my $sub ( sort keys %subref ) 
	{
		my $refbase = substr $sub, 0, 1;		
		my $mutantbase = substr $sub, 2, 1;

		# skip if middle base of trinuc is not same as current sub reference base
		if ( $refbase ne $trirefbase )
		{
			next;
		}

		my $lag = 0;
		my $lead = 0;
		my $int = 0;
		if ( exists $lagbin{$tri}{$sub} )
		{
			$lag = $lagbin{$tri}{$sub};
		}
		if ( exists $leadbin{$tri}{$sub} )
		{
			$lead = $leadbin{$tri}{$sub};
		}
                if ( exists $intergenicbin{$tri}{$sub} )
                {
                        $int = $intergenicbin{$tri}{$sub};
                }
                print OUT "$tri\t$mutantbase\t$sub\t$lag\t$lead\t$int\n";
	}
}

