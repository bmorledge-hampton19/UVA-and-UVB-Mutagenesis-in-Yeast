#!/usr/bin/perl

use strict;
use warnings;

package OriginCoord;

sub new
{
	my ($class) = @_;
	
	my $self = bless {}, $class;

	# open file with gene positions
	open( ORIGIN, "../Origins.txt") || die "Couldn't open file\n";
	my $header = <ORIGIN>;
	
	my %chromosome;
	my %originstart;
	my %originend;
	my %strand;
	
=pod
	my %chr_convert = ( "chrI" => "chr1","chrII" => "chr2", "chrIII" => "chr3", "chrIV" => "chr4","chrV" => "chr5", "chrVI" => "chr6", "chrVII" => "chr7", "chrVIII" => "chr8", "chrIX" => "chr9", "chrX" => "chr10", "chrXI" => "chr11", "chrXII" => "chr12", "chrXIII" => "chr13", "chrXIV" => "chr14", "chrXV" => "chr15", "chrXVI" => "chr16" ); 
=cut

	while( my $line = <ORIGIN> )
	{
		chomp($line);
		my @fields = split /\t/, $line;
		if( $fields[1] =~ /^(ARS[0-9]+\.*[0-9]*)/ )
		{
			my $acc = $1;
			$originstart{$acc} = $fields[2]; 
			$originend{$acc} = $fields[3];
			push @{$chromosome{$fields[0]}}, $acc ;
		}
	}

	close ( ORIGIN );
	$self->{'chromosome'} = \%chromosome;
	$self->{'originstart'} = \%originstart;
	$self->{'originend'} = \%originend;

	return $self;	
}

sub get_chromosomes
{
	my ($self) = @_;

	return %{$self->{'chromosome'}};

}

sub get_origin_starts
{
        my ($self) = @_;

        return %{$self->{'originstart'}};
}

sub get_origin_ends
{
        my ($self) = @_;

        return %{$self->{'originend'}};
}

	
	

1;	
