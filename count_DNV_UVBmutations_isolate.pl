#!/usr/bin/perl -w

use strict;
use warnings;

my $header = <STDIN>;
print STDERR "Counting isolates - Header: $header\n";

my %isolate = ();
my %isolatekey = ();
my %refbasekey = ();
my %subkey = ();
my $mitocount = 0;

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


=pod
	#add run number
	if ( $fields[5] =~ /^[12]$/ )
	{
		$genotype .= "_run$fields[5]";
	}
	else
	{
		die "Error with line for run: $line\n";
	}

	# remove Rad16_25 and Rad26_17 isolate (low sequencing depth)
	if ( $fields[0] eq "Rad26_17" || $fields[0] eq "dRad16_25" )
	{
		print STDERR "skipped line: $line\n";
		next;
	}
=cut

        # remove B01 Rad16 UVB and C13 Rad26 UVB isolate (presumably low sequencing depth)
        if ( $fields[0] eq "B01 Rad16 UVB" || $fields[0] eq "C13 Rad26 UVB" )
        {
                print STDERR "skipped line: $fields[0]\n";
                $skippedlines++;
                next;
        }

        my $chr = $fields[2];

        if ( $chr eq "chrMito" )
        {
                print STDERR "Skipped line: ${chr}\n";
                $skippedlines++;
                next;
        }

	if ( $fields[11] =~ /MNV/ )
	{
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
		elsif ( length( $origbase ) != 2 || length( $newbase ) != 2 )
		{
			print STDERR "non-2bp MNV: $line\n";
			next;
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
		$refbasekey{$refbase} = 1;
		my $sub = $refbase . ">" . $mutbase;
		$subkey{$sub} = 1;
		$isolatekey{$genotype}{$fields[0]}++;
		$isolate{$genotype}{$refbase}{$sub}{$fields[0]}++;

	}
}

my $strain = "";
my $basename = "_isolate_DNVcount.txt";
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
			print OUT "\tReference base\tMutant base\tMutation substitution\tMean\tSEM\tN\tTotal Mutations\tStandard deviation";
			print OUT "\n";

                }
        }
	foreach my $refbase ( sort keys %refbasekey )
	{
	        foreach my $sub ( sort keys %subkey )
  	      	{
			my $ref = substr $sub, 0, 2;
			if ( $ref ne $refbase )
			{
				next;
			}

                        my $middle = substr $sub, 2, 1;
			if ( length($sub) != 5 || $middle ne ">" )
			{
				die "Error with sub: $sub\n";
			}
			my $mutbase = substr $sub, 3, 2;

			if ( $refbase eq $mutbase )
			{
				next;
			}
			print OUT "$refbase\t$mutbase\t$sub";

			my @mutvals = ();

			foreach my $iso ( sort keys %{$isolatekey{$strain}} )
			{
				my $mutcount = 0;
				if ( exists $isolate{$gene}{$refbase}{$sub}{$iso} )
				{
					$mutcount = $isolate{$gene}{$refbase}{$sub}{$iso};
				}
				print OUT "\t$mutcount";
				push @mutvals, $mutcount;

			}
			my $sum = &summation(\@mutvals);
			my $avg = &average(\@mutvals);
			my $stdev = &stdev(\@mutvals);
			my $sem = &stderror(\@mutvals);
			my $n = scalar @mutvals;
			print OUT "\t$refbase\t$mutbase\t$sub\t$avg\t$sem\t$n\t$sum\t$stdev";
			print OUT "\n";
		}
	}
	
}

close ( OUT );

print STDERR "Mitochondrial mutations excluded: $mitocount\n";

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
