#!/usr/bin/perl -w

use strict;
use warnings;

my $header = <STDIN>;
print STDERR "SNV processing\tHeader: $header\n";
my %isolate = ();
my %isolatekey = ();
my %refkey = ();
my %mutkey = ();
my %subkey = ();
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

        if ( $id =~ /([a-zA-Z0-9]+) ([A-Za-z]+) ([NOUV]+)/ )
        {
                $isolatename = $1;
                $genotype = $2;

                if ( $genotype eq "OGG" )
                {
                        $genotype = "ogg1";
                }
                my $uv = $3;
                if ( $uv eq "NO" )
                {
                        $uv = "noUVA";
                }
                elsif ( $uv eq "UV" )
                {
                        $uv = "UVA";
                }
                else
                {
                        die "problem with name\n";
                }
                $uv = "_" . $uv;
                $genotype .= $uv;
        }
        else
        {
                die "error with id: $id\n";
        }

        # remove D15 Ogg1d isolate (low sequencing depth)
        if ( $fields[0] eq "D15 OGG UV" )
        {
                print STDERR "skipped line: $fields[0]\n";
                $skippedlines++;
                next;
        }

	$isolatekey{$genotype}{$fields[0]}++;

	if ( $fields[9] =~ /SNV/ )
	{
        	my $chr = $fields[2];

                if ( $chr eq "chrMito" )
                {
                        print STDERR "Skipped line: ${chr}\n";
			$skippedlines++;
			next;
                }

	        my $mut_end = $fields[3];

        	# convert to 0-based bedfile format from 1-based mutation call
	        my $mut_start = $mut_end - 1;

		# Get strand and base information for mutation
		my $origbase = $fields[11];
		my $newbase = $fields[12]; 
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

		$mutkey{$mutbase} = 1;
		$refkey{$refbase} = 1;
		$isolate{$genotype}{$refbase}{$mutbase}{$fields[0]}++;
	}
}
my $strain = "";
my $basename = "_mutype_isolate_Mutcount.txt";
my $filename = "";

foreach my $gene ( sort keys %isolate )
{
        if ( $gene ne $strain )
        {
                if ( $gene =~ /[A-Za-z0-9]+/ )
                {
                        if ( $strain ne "" )
                        {
				close ( OUT );
                        }
                        $strain = $gene;
                        $filename = $strain . $basename;
                        open ( OUT, ">$filename" ) or die "cannot open file: $filename\n";
			print OUT "Isolate Mutation count for $gene\n";
			print OUT "Reference base\tMutant base\tMutation substitution";
			my $num = 1;
			foreach my $iso ( sort keys %{$isolatekey{$strain}} )
			{
				print OUT "\t${iso} (#$num)";
				$num++;
			} 	
			print OUT "\tReference base\tMutation substitution\tMean\tSEM\tN\tTotal Mutations\tStandard deviation";
			print OUT "\n";

                }
        }
	foreach my $refbase ( sort keys %refkey )
	{
	        foreach my $mutbase ( sort keys %mutkey )
  	      	{
			if ( $refbase eq $mutbase )
			{
				next;
			}
                	my $sub = $refbase . ">" . $mutbase;

			$subkey{$sub} = 1;
			print OUT "$refbase\t$mutbase\t$sub";

			my @mutvals = ();

			foreach my $iso ( sort keys %{$isolatekey{$strain}} )
			{
				my $mutcount = 0;
				if ( exists $isolate{$gene}{$refbase}{$mutbase}{$iso} )
				{
					$mutcount = $isolate{$gene}{$refbase}{$mutbase}{$iso};
				}
				print OUT "\t$mutcount";
				push @mutvals, $mutcount;

			}
			my $sum = &summation(\@mutvals);
			my $avg = &average(\@mutvals);
			my $stdev = &stdev(\@mutvals);
			my $sem = &stderror(\@mutvals);
			my $n = scalar @mutvals;
			print OUT "\t$refbase\t$sub\t$avg\t$sem\t$n\t$sum\t$stdev";
			print OUT "\n";
		}
	}
	
}

close ( OUT );

print STDERR "Skipped lines: $skippedlines\n";

# code adapted from https://edwards.sdsu.edu/research/calculating-the-average-and-standard-deviation/
sub average
{
        my($data) = @_;
        if (not @$data) {
                die("Empty arrayn");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = 1.0 * $total / @$data;
        return $average;
}
sub stdev
{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}
sub stderror
{
        my($data) = @_;
	my $stdval = &stdev($data);

	my $semval = $stdval / ( ( @$data ) ** 0.5 );
        return $semval;
}
sub summation
{
        my($data) = @_;
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
	return $total;
} 
