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

my $chr = "";
my %strandlookup;
my $mitomutcount = 0;

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
		my $leftend = $oristart{$acc} - 1;
		# strand indicates leading strand
		print "$chr\t$leftstart\t$leftend\t$acc\t.\t$leftstrand\n";

		# indicate strand of leading strand --> right of origin, - strand is leading
		my $rightstrand = "-";
		my $rightstart = $oriend{$acc}; # because of bed formatting
		my $rightend = $oriend{$acc} + $window;
		print "$chr\t$rightstart\t$rightend\t$acc\t.\t$rightstrand\n";

	}

	print STDERR "$skipped skipped origins\n";
}
