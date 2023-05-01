#!/usr/bin/perl

use strict;
use warnings;

package GeneCoord;

sub new
{
	my ($class) = @_;
	
	my $self = bless {}, $class;

	# open file with gene positions
	open( GENE, "../TSS_Iyer_NAR2014_saccer3_1based.txt") || die "Couldn't open file\n";
	my $header = <GENE>;
	
	my %chromosome;
	my %tss;
	my %tts;
	my %strand;

	# Exclude genes overlapping with rDNA
	my @excluded = ("YLR154C-G", "YLR154C-H", "YLR154W-C", "YLR155C", "YLR161W", "YLR162W", "YLR162W-A");	

	# for Rad16 CPD-seq experiments, exclude Rad16 
	#push @excluded, "YBR114W";

        # for Rad26 CPD-seq experiments, exclude Rad16
        #push @excluded, "YJR035W";

	# exclude cup1-1 related genes, seems to be amplified
	push @excluded, "YHR055C";
        push @excluded, "YHR054C";
        push @excluded, "YHR053C";

	my %exclude_acc = ();
	foreach my $exc (@excluded)
	{
		$exclude_acc{$exc} = 1;
	}

	# Polyadenylation site date (e.g., tts)
        open (PAS, "../PAS_Iyer_NAR2014_saccer3_1based.txt" ) || die "couldn't open file\n";
	$header = <PAS>;

	while( my $line = <GENE> )
	{
		chomp($line);
		my @fields = split /\t/, $line;
		if( $fields[2] =~ /^(Y[A-P][LR][0-9]{3}[CW]\-?[A-H]?)/ )
		{
			my $acc = $1;
			# remove rDNA overlapping genes
			if ( $exclude_acc{$acc} )
			{
				print STDERR "$acc is excluded\n";
				next;
			}
			$tss{$acc} = $fields[1]; 
			push @{$chromosome{$fields[0]}}, $acc;

			my $str = "";
			if ( $acc =~ /^Y[A-P][LR][0-9]{3}W/ )
			{
				$str = "+";
			}
			elsif ( $acc =~ /^Y[A-P][LR][0-9]{3}C/ )
			{
				$str = "-";
			}
			else
			{
				die "No strand information for gene $acc\n";
			}

			$strand{$acc} = $str;
		}
		else
		{
			print STDERR "accession doesn't fit pattern: $fields[2]\n";
		}
	}

	close ( GENE );

        while( my $line = <PAS> )
        {
                chomp($line);
                my @fields = split /\t/, $line;
                if( $fields[2] =~ /^(Y[A-P][LR][0-9]{3}[CW]\-?[A-H]?)/ )
                {
                        my $acc = $1;
                        # remove rDNA overlapping genes
                        if ( $exclude_acc{$acc} )
                        {
                                print STDERR "$acc is excluded\n";
                                next;
                        }
                        if ( exists $tss{$acc} )
			{
				$tts{$acc} = $fields[1];
			} 
                        else
                        {
                                #print STDERR "No tss information for gene $acc\n";
                        }
                }
        }

	close ( PAS );

	$self->{'chromosome'} = \%chromosome;
	$self->{'tss'} = \%tss;
	$self->{'tts'} = \%tts;
	$self->{'strand'} = \%strand;

	return $self;	
}

sub get_chromosomes
{
	my ($self) = @_;

	return %{$self->{'chromosome'}};

}

sub get_tss
{
        my ($self) = @_;

        return %{$self->{'tss'}};
}

sub get_tts
{
        my ($self) = @_;

        return %{$self->{'tts'}};
}

sub get_strand
{
        my ($self) = @_;

        return %{$self->{'strand'}};
}
	
	

1;	
