#!/usr/bin/perl -w
## Author: Ashok K. Sharma
use strict;


#Thioalkalivibrio_nitratireducens        k_Bacteria;p_Proteobacteria;f_Ectothiorhodospiraceae;g_Thioalkalivibrio;s_sulfidophilus
#Thermodesulfobium_narugense     k_Bacteria;p_Firmicutes;f_Thermodesulfobiaceae;g_Thermodesulfobium;s_narugense

my $input = <$ARGV[0]>;
open (IN, $input) || die "Cannot open the input file\n";

my $hash_taxa={};
my $hash={};


while (chomp(my $line = <IN>))
{
	my @x = split("\t", $line);
	my $Species = $x[0];
	my $Taxa = $x[1];
	my @T = split (";", $Taxa);
	(my $genus) = $line =~ /;g_(.*);s/;
	(my $family) = $line =~ /;f_(.*);g/;
	(my $phylum) = $line =~ /;p_(.*);f/;
	$hash_taxa->{$Species}->{'GENUS'}=$genus;
	$hash_taxa->{$Species}->{'FAMILY'}=$family;
	$hash_taxa->{$Species}->{'PHYLUM'}=$phylum;
}	

my $input2 = <$ARGV[1]>;
open (IN2, $input2) || die "Cannot open the input2 file\n";

# 411489.CLOL250_01016	Clostridium_sp.	100.00
# 411489.CLOL250_01016	Clostridium_sp.	100.00

while (chomp (my $line2 = <IN2>))
{
	my @y = split("\t", $line2);
	push (@{$hash->{$y[0]}->{'Target'}},$y[1]);
	push (@{$hash->{$y[0]}->{'identity'}},$y[2]);
}

open (OUT, ">Taxonomy_Assignment_LCA_Hits") || die "Cannot open the output file\n";

open (OUT2, ">Taxonomy_Unannotated_LCA_Hits") || die "Cannot open the file\n";

foreach my $k (keys %{$hash})
{
	my $M = &Comp(@{$hash->{$k}->{'Target'}});
	if ($M == 0)
	{
		print OUT $k."\t".$hash->{$k}->{'Target'}->[0]."\t".$hash->{$k}->{'identity'}->[0]."\t"."S"."\n";
#		print "yes Species are equal";<stdin>;
	}
	else 
	{
		my @Genus;
		foreach my $s (@{$hash->{$k}->{'Target'}})
		{
			if (exists($hash_taxa->{$s}->{'GENUS'}))
			{
				push (@Genus, $hash_taxa->{$s}->{'GENUS'});
#				print "Genus exists";<stdin>;
			}
			else 
			{
				print OUT2 $k."\t".$s."\n";
			}
		}

		if (scalar(@Genus) > 0)
		{
			$M = &Comp (@Genus);
		}
		else 
		{
			next
		}

		if ($M == 0)
		{
			print OUT $k."\t".$Genus[0]."\t".$hash->{$k}->{'identity'}->[0]."\t"."G"."\n";
#				print "Genus are equal";<stdin>;
		}
		else 
		{
			my @Family;	
			foreach my $g (@{$hash->{$k}->{'Target'}})
			{
				if (exists($hash_taxa->{$g}->{'FAMILY'}))
				{
					push (@Family, $hash_taxa->{$g}->{'FAMILY'});
#						print $hash_taxa->{$g}->{'FAMILY'};<stdin>;
				}
				else 
				{
					print OUT2 $k."\t".$g."\n";
				}
			}

			if (scalar(@Family) > 0)
			{	
				$M = &Comp (@Family);
			}
			else 
			{
				next
			}

			if ($M == 0)
			{
				print OUT $k."\t".$Family[0]."\t".$hash->{$k}->{'identity'}->[0]."\t"."F"."\n";
			}
			else 
			{
				my @Phylum;
				foreach my $f (@{$hash->{$k}->{'Target'}})
				{
					if (exists($hash_taxa->{$f}->{'PHYLUM'}))
					{
#								print "Phylum exists";
						push (@Phylum, $hash_taxa->{$f}->{'PHYLUM'});
#								print $hash_taxa->{$f}->{'PHYLUM'};<stdin>;
					}
					else {
						print OUT2 $k."\t".$f."\n";
						
						}
				}
				if (scalar(@Phylum) > 0)
				{
					$M = &Comp(@Phylum);
				}
				else 
				{
					next
				}
				
				if ($M == 0)
				{
						print OUT $k."\t".$Phylum[0]."\t".$hash->{$k}->{'identity'}->[0]."\t"."P"."\n";
				}
			}
		}
	}
}


	


sub Comp {
		my ($taxa) = shift (@_);
		my $N = 0;
		foreach (@_)
		{
			if ($_ ne $taxa)
			{
				$N++
			}
		}
		$N;
}
