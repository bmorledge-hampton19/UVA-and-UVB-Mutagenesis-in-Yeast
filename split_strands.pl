#/usr/bin/perl

use strict;
use warnings;

print "Input name of bed file:\n";
my $filename = <STDIN>;
chomp $filename;
open( FILE, "$filename" ) || die "Couldn't open file $filename\n";

$filename =~ s/\.bed//;
 
my $plusfile = $filename . "_plusstrand.bed";
my $minusfile = $filename . "_minusstrand.bed";

open( PLUS, ">$plusfile" );
open( MINUS, ">$minusfile" );
while ( my $line = <FILE> )
{
	if( $line =~ /\+/ )
	{
		print PLUS $line;
	}
	elsif( $line =~ /-/ )
	{
		print MINUS $line;
	}
	else 
	{	 
		print "Error: line $line is neither plus or minus\n";	
	}
}
