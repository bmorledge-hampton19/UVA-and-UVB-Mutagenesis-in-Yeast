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
my %isolate = ();
my %isolatekey = ();
my %trinuckey = ();
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

		$mutkey{$mutbase} = 1;
		$trinuckey{$trinuc} = 1;
		$isolatekey{$genotype}{$fields[0]}++;
		$isolate{$genotype}{$trinuc}{$mutbase}{$fields[0]}++;
	}
}
my $strain = "";
my $basename = "_trinuc_isolate_Mutcount.txt";
my $filename = "";
my $classname = "_trinuclass_isolate_mutcount.txt";

my %triclass = ();
my %classkey = ( "Dipy" => "Dipyrimidine", "five" => "5\'", "three" => "3\'", "z No" => "No");

foreach my $gene ( sort keys %isolate )
{
        if ( $gene ne $strain )
        {
                if ( $gene =~ /[A-Za-z0-9]+/ )
                {
                        if ( $strain ne "" )
                        {
				close ( OUT );
				# open new file with class counts
	                        my $classfile = $strain . $classname;
        	                open ( CLASS, ">$classfile" ) or die "cannot open file: $classfile\n";
                	        print CLASS "Isolate Mutation class count for $gene\n";
                        	print CLASS "Mutation substitution\tClass";
	                        my $num = 1;
        	                foreach my $iso ( sort keys %{$isolatekey{$strain}} )
                	        {
                        	        print CLASS "\t${iso} (#$num)";
                                	$num++;
	                        }
        	                print CLASS "\tMutation substitution\tClass\tMean\tSEM\tN\tTotal Mutations\tStandard deviation";
                	        print CLASS "\n";
			        foreach my $sub ( sort keys %subkey )
			        {
		                        print CLASS "$sub";

					foreach my $c ( sort keys %classkey )
					{						
						print CLASS "\t$classkey{$c}";
						my @mutvals = ();
						foreach my $iso ( sort keys %{$isolatekey{$strain}} )
						{
							my $mutcount = 0;
							if ( exists $triclass{$sub}{$c}{$iso} )
							{
								$mutcount = $triclass{$sub}{$c}{$iso};
							}	
							print CLASS "\t$mutcount";
							push @mutvals, $mutcount;
		                                }
		                        	my $sum = &summation(\@mutvals);
			                        my $avg = &average(\@mutvals);
			                        my $stdev = &stdev(\@mutvals);
			                        my $sem = &stderror(\@mutvals);
			                        my $n = scalar @mutvals;
		        	                print CLASS "\t$sub\t$classkey{$c}\t$avg\t$sem\t$n\t$sum\t$stdev";
		                	        print CLASS "\n";
					}
				}
				close ( CLASS );

				#reset triclass variable
				%triclass = ();

                        }
                        $strain = $gene;
                        $filename = $strain . $basename;
                        open ( OUT, ">$filename" ) or die "cannot open file: $filename\n";
			print OUT "Isolate Mutation count for $gene\n";
			print OUT "Trinucleotide context\tMutant base\tMutation substitution";
			my $num = 1;
			foreach my $iso ( sort keys %{$isolatekey{$strain}} )
			{
				print OUT "\t${iso} (#$num)";
				$num++;
			} 	
			print OUT "\tTrinucleotide context\tMutation substitution\tMean\tSEM\tN\tTotal Mutations\tStandard deviation";
			print OUT "\n";

                }
        }
	foreach my $tri ( sort keys %trinuckey )
	{
	        foreach my $mutbase ( sort keys %mutkey )
  	      	{
               		my $refbase = substr $tri, 1, 1;
			if ( $refbase eq $mutbase )
			{
				next;
			}
                	my $sub = $refbase . ">" . $mutbase;

			$subkey{$sub} = 1;
			print OUT "$tri\t$mutbase\t$sub";

			my @mutvals = ();

	                # find trinuc class (e.g., 5' dipy, etc.)
                        my $class = "";
			my $dipy = 0;

                        if ( $tri =~ /^[CT][CT][CT]$/ )
			{
				$class = "both";
				$dipy = 1;
			}
			elsif ( $tri =~ /^[AG][CT][CT]$/ )
			{
				$class = "five";
				$dipy = 1;
			}
			elsif ( $tri =~ /^[CT][CT][AG]$/ )
			{
				$class = "three";
				$dipy = 1;
			}
			elsif ( $tri =~ /^[AG][CT][AG]$/ )
			{
				$class = "z No";
				$dipy = 0;
			}
			else
			{
				die "Error with trinuc: $tri\n";
			}

			foreach my $iso ( sort keys %{$isolatekey{$strain}} )
			{
				my $mutcount = 0;
				if ( exists $isolate{$gene}{$tri}{$mutbase}{$iso} )
				{
					$mutcount = $isolate{$gene}{$tri}{$mutbase}{$iso};
				}
				print OUT "\t$mutcount";
				push @mutvals, $mutcount;

				if ( $dipy )
				{
					$triclass{$sub}{"Dipy"}{$iso} += $mutcount;
				}	

				if ( $class eq "both" )
				{
					$triclass{$sub}{"five"}{$iso} += $mutcount;
					$triclass{$sub}{"three"}{$iso} += $mutcount;
				}
				else
				{
					$triclass{$sub}{$class}{$iso} += $mutcount;
				}
			}
			my $sum = &summation(\@mutvals);
			my $avg = &average(\@mutvals);
			my $stdev = &stdev(\@mutvals);
			my $sem = &stderror(\@mutvals);
			my $n = scalar @mutvals;
			print OUT "\t$tri\t$sub\t$avg\t$sem\t$n\t$sum\t$stdev";
			print OUT "\n";
		}
	}
	
}

close ( OUT );
# open new file with class counts
my $classfile = $strain . $classname;
open ( CLASS, ">$classfile" ) or die "cannot open file: $classfile\n";
print CLASS "Isolate Mutation class count for $strain\n";
print CLASS "Mutation substitution\tClass";
my $num = 1;
foreach my $iso ( sort keys %{$isolatekey{$strain}} )
{
        print CLASS "\t${iso} (#$num)";
        $num++;
}
print CLASS "\tMutation substitution\tClass\tMean\tSEM\tN\tTotal Mutations\tStandard deviation";
print CLASS "\n";
foreach my $sub ( sort keys %subkey )
{
        print CLASS "$sub";

        foreach my $c ( sort keys %classkey )
        {
                print CLASS "\t$classkey{$c}";
                my @mutvals = ();
                foreach my $iso ( sort keys %{$isolatekey{$strain}} )
                {
                        my $mutcount = 0;
                        if ( exists $triclass{$sub}{$c}{$iso} )
                        {
                                $mutcount = $triclass{$sub}{$c}{$iso};
                        }
                        print CLASS "\t$mutcount";
                        push @mutvals, $mutcount;
                }
                my $sum = &summation(\@mutvals);
                my $avg = &average(\@mutvals);
                my $stdev = &stdev(\@mutvals);
                my $sem = &stderror(\@mutvals);
                my $n = scalar @mutvals;
                print CLASS "\t$sub\t$classkey{$c}\t$avg\t$sem\t$n\t$sum\t$stdev";
                print CLASS "\n";
        }
}
close ( CLASS );

print STDERR "Skipped lines: $skippedlines\n";
print STDERR "Mismatched trinucs: $mismatchtri\n";

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
