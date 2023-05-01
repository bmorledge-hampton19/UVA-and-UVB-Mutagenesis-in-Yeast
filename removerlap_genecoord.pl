#!/usr/bin/perl

use strict;
use warnings;

my @prevline = ();

while ( my $line = <STDIN> )
{
	chomp $line;
	my @currline = split /\t/, $line;

	# check for merge: overlapping genes - assume file is sorted based on chr and start and end pos
	# note: not <= because this is a bed file, end position is 1-based, not 0-based
	
	if( $currline[1] < $prevline[2] && $currline[0] eq $prevline[0] )
	{
		# make sure on different strands -- should be after previous merge
		if ( $currline[5] eq $prevline[5] )
		{
			die "error same strand!!!\n";
		}
		elsif ( $prevline[2] > $currline[2] )
		{
			die "Complete overlap prevline: @prevline\tCurrent line: @currline\n";
		}
		else
		{
			my $temp = $currline[1];
			$currline[1] = $prevline[2];
			$prevline[2] = $temp;

			# update coordinates in header
			my $prevstart = $prevline[1] + 1;
			$prevline[3] = "$prevline[0]:$prevstart-$prevline[2]";

			my $currstart = $currline[1] + 1;
			$currline[3] = "$currline[0]:$currstart-$currline[2]";

              	  	if ( exists $prevline[0] )
                	{
                        	my $last = (scalar @prevline ) - 1;
	                        # print prev line
        	                for ( my $i = 0; $i < $last; $i++ )
                	        {
                        	        print "$prevline[$i]\t";
                        	}
                        	# print last field
                        	print "$prevline[$last]\n";
                	}
                	@prevline = @currline;
		}
	}
	else
	{
		if ( exists $prevline[0] )
		{
                	my $last = (scalar @prevline ) - 1;
			# print prev line
			for ( my $i = 0; $i < $last; $i++ )
			{
				print "$prevline[$i]\t";
			}
			# print last field 
			print "$prevline[$last]\n";
		}

		@prevline = @currline;
	}

}

# print line remaining in prevline
my $last = (scalar @prevline ) - 1;
# print prev line 
for ( my $i = 0; $i < $last; $i++ )
{
        print "$prevline[$i]\t";
}
# print last field 
print "$prevline[$last]\n";

