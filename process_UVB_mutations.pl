#!/usr/bin/perl -w

use strict;
use warnings;

my $fafile = "trinucontext_UVB_mutations_yeast.fa";
my %tricontext = ();

# process trinuc values
open ( FA, $fafile ) || die "couldn't open $fafile\n";
while ( my $line1 = <FA> )
{       
        chomp ( $line1 );
        $line1 =~ s/^>// || die "misformatted line: $line1\n";
        my $trival = <FA> || die "missing line\n";
        chomp ( $trival );
        $tricontext{$line1} = $trival;
}

my $header = <STDIN>;
print STDERR "SNV processing\tHeader: $header\n";
my $strain = "";
my $basename = "_Muts_allSNVs.bed";
my $filename = "";
my $filteredcount = 0;
my $mismatchtri = 0;
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
				print STDERR "$strain: Mismatched trinucs = $mismatchtri\n";
				$mismatchtri = 0;
				close ( OUT );
			}
			$strain = $genotype;
			$filename = $strain . $basename;
			open ( OUT, ">$filename" ) or die "cannot open file: $filename\n";
		}	
	}			
	if ( $fields[11] =~ /SNV/ )
	{
        	my $chr = $fields[2];

                if ( $chr eq "chrMito" )
                {
                        print STDERR "Skipped line: ${chr}\n";
			$skippedlines++;
			next;
                }

	        my $mut_end = $fields[4];

        	# convert to 0-based bedfile format from 1-based mutation call
	        my $mut_start = $mut_end - 1;

		if ( $mut_end != $fields[3] )
		{
			die "mismatched start and end position!";
		}

		# Get strand and base information for mutation
		my $origbase = $fields[12];
		my $newbase = $fields[13]; 
		my $strand = "";
		my $refbase = "";
		my $mutbase = "";
		my $trinuc = "";
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

		# fix trinuc
		my $tristart = $mut_start - 1;
		my $triend = $mut_end + 1;
		my $trikey = "${chr}:$tristart-${triend}(${strand})";

		if ( exists $tricontext{$trikey} )
		{
			$trinuc = $tricontext{$trikey};
		}
		else
		{
			die "no trinuc value for $trikey\n\n";
		}

		my $midtri = substr $trinuc, 1, 1;

		if ( $refbase ne $midtri )
		{
			die "mismatch in mutated base in reference: $refbase\t$midtri\n";
		}

        	print OUT "${chr}\t$mut_start\t$mut_end\t$trinuc\t$mutbase\t$strand\n";

	}
	else 
	{
		$filteredcount++;
		#print STDERR "Filtered line: $line\n";
	}
}
print STDERR "Filtered lines: $filteredcount\n";
print STDERR "Skipped lines: $skippedlines\n\n";
