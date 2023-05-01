#!/usr/bin/perl -w

use strict;
use warnings;

my $header = <STDIN>;
print STDERR "DNP mutation processing - Header: $header\n";
my $strain = "";
my $basename = "_Muts_allDNVs.bed";
my $filename = "";
my $filteredcount = 0;
my $skippedlines = 0;
while( my $line = <STDIN> )
{
	chomp $line;
	if ( $line eq "")
	{
		next;
	}
	my @fields = split /\t/, $line;

        my $id = $fields[0];
        my $genotype = "";
        my $isolatename = "";
        if ( $id =~ /([a-zA-Z0-9]+) ([A-Za-z0-9]+\s?R?P?) UVB/ )
        {
                $isolatename = $1;
                $genotype = $2;
                $genotype .= "_UVB";
        }
        else
        {
                die "error with id: $id\n";
        }

        # remove B01 Rad16 UVB and C13 Rad26 UVB isolate (presumably low sequencing depth)
        if ( $fields[0] eq "B01 Rad16 UVB" || $fields[0] eq "C13 Rad26 UVB" )
        {
                print STDERR "skipped line: $fields[0]\n";
                $skippedlines++;
                next;
        }

	if ( $genotype ne $strain )
	{
		if ( $genotype =~ /[A-Za-z0-9]+/ )
		{
			if ( $strain ne "" )
			{
				close ( OUT );
			}
			$strain = $genotype;
			$filename = $strain . $basename;
			open ( OUT, ">$filename" ) or die "cannot open file: $filename\n";
		}	
	}			
	if ( $fields[11] =~ /MNV/ )
	{
        	my $chr = $fields[2];

                if ( $chr eq "chrMito" )
                {
                        print STDERR "Skipped line: Chr${chr}\n";
                        $skippedlines++;
                        next;
                }

	        my $mut_end = $fields[4];

        	# convert to 0-based bedfile format from 1-based mutation call
	        my $mut_start = $fields[3] - 1;

		if ( $mut_end - $mut_start != 2 )
		{
                        print STDERR "Skipped line: Chr${chr}\n";
                        $skippedlines++;
                        next;
		}
		# Get strand and base information for mutation
		my $origbase = $fields[12];
		my $tri = $origbase;
		my $newbase = $fields[13]; 
		my $strand = "";
		my $refbase = "";
		my $mutbase = "";
		if ( $origbase eq $newbase )
		{
			die "Reference base equals mutbase: $line\n";
		}
		if ( $origbase =~ /[CT][CT]/ || $origbase =~ /AC/ || $origbase =~ /CA/ )
		{
			$strand = "+";
			$refbase = $origbase;
			$mutbase = $newbase;
		}
		elsif ( $origbase =~ /[AG][AG]/ || $origbase =~ /GT/ || $origbase =~ /TG/ )
		{
			$strand = "-";
			$refbase = reverse $origbase;
			$refbase =~ tr/ATGC/TACG/;
			$mutbase = reverse $newbase;
			$mutbase =~ tr/ATGC/TACG/;
		}
		elsif ( $origbase =~ /AT/ || $origbase =~ /TA/ || $origbase =~ /CG/ || $origbase =~ /GC/ )
		{
			# use mutbase to stratify
	                if ( $newbase =~ /[CT][CT]/ || $newbase =~ /AC/ || $newbase =~ /CA/ )
	                {
	                        $strand = "+";
	                        $refbase = $origbase;
	                        $mutbase = $newbase;
	                }
	                elsif ( $newbase =~ /[AG][AG]/ || $newbase =~ /GT/ || $newbase =~ /TG/ )
	                {
	                        $strand = "-";
	                        $refbase = reverse $origbase;
	                        $refbase =~ tr/ATGC/TACG/;
	                        $mutbase = reverse $newbase;
	                        $mutbase =~ tr/ATGC/TACG/;
	                }
			else
			{	
				# assign to plus strand by default
				$strand = "+";
                        	$refbase = $origbase;
                        	$mutbase = $newbase;
			}
		}
		else
		{
			die "Error if double mutant: $origbase\n";
		}
        	print OUT "chr${chr}\t$mut_start\t$mut_end\t$refbase\t$mutbase\t$strand\n";
	}
	else 
	{
		$filteredcount++;
		#print STDERR "Filtered line: $line\n";
	}
}
print STDERR "Filtered lines: $filteredcount\n";
print STDERR "SKIPPED LINES: $skippedlines\n\n";
