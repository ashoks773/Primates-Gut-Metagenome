#!/usr/bin/perl
use strict;


#Thioalkalivibrio_nitratireducens        k_Bacteria;p_Proteobacteria;f_Ectothiorhodospiraceae;g_Thioalkalivibrio;s_sulfidophilus
#Thermodesulfobium_narugense     k_Bacteria;p_Firmicutes;f_Thermodesulfobiaceae;g_Thermodesulfobium;s_narugense

my $input = <$ARGV[0]>;
open (IN, $input) || die "Cannot open the input file\n";

my $hash_taxa={};
my $hash={};

while (chomp (my $line = <IN>))
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
	$hash_taxa->{$genus}->{'GENPHY'}=$phylum;
	$hash_taxa->{$family}->{'FAMPHY'}=$phylum;
}	

my $input2 = <$ARGV[1]>;
open (IN2, $input2) || die "Cannot open the input2 file\n";

# 411489.CLOL250_01016	Clostridium_sp.	100.00	S
# SZEY-101A_GL0027349	Faecalibacterium_cf.	99.21	S
# Gene_ID_UA_gene_586439	Oscillibacter_valericigenes	82.50	S

while (chomp (my $line2 = <IN2>))
{
	my @y = split("\t", $line2);
	$hash->{$y[0]}->{'Target'}=$y[1];
	$hash->{$y[0]}->{'identity'} = $y[2];
	$hash->{$y[0]}->{'Ass'}=$y[3];
}

open (OUT, ">GENES_FINAL_TAXA_ASSIGNMENT") || die "Cannot open the output file\n";

open (OUT2, ">GENES_FINAL_TAXA_UNASSIGNED") || die "Cannot open the file\n";

foreach my $k (keys %{$hash})
{
	if ($hash->{$k}->{'identity'} >= 95)
	{
		print OUT $k."\t".$hash->{$k}->{'Target'}."\n";
	}
	elsif (($hash->{$k}->{'identity'} < 95) && ($hash->{$k}->{'identity'} >= 85))
	{
		
		if ($hash->{$k}->{'Ass'} eq "G")
		{
			print OUT $k."\t".$hash->{$k}->{'Target'}."\n";
		}
		elsif ($hash->{$k}->{'Ass'} eq "F")
		{
			print OUT $k."\t".$hash->{$k}->{'Target'}."\n";
		}
		elsif ($hash->{$k}->{'Ass'} eq "P")
		{
			print OUT $k."\t".$hash->{$k}->{'Target'}."\n";
		}
		elsif ($hash->{$k}->{'Ass'} eq "S")
		{
			if (exists($hash_taxa->{$hash->{$k}->{'Target'}}->{'GENUS'}))
			{
				print OUT $k."\t".$hash_taxa->{$hash->{$k}->{'Target'}}->{'GENUS'}."\n";
			}
		}  
		else 
		{
			print OUT2 $k."\t".$hash->{$k}->{'Target'}."\n";
		}
	}
	elsif (($hash->{$k}->{'identity'} < 85) && ($hash->{$k}->{'identity'} >=65))
	{
		if ($hash->{$k}->{'Ass'} eq "G")
                {
                        print OUT $k."\t".$hash_taxa->{$hash->{$k}->{'Target'}}->{'GENPHY'}."\n";
                }
                elsif ($hash->{$k}->{'Ass'} eq "F")
                {
                        print OUT $k."\t".$hash_taxa->{$hash->{$k}->{'Target'}}->{'FAMPHY'}."\n";
                }
                elsif ($hash->{$k}->{'Ass'} eq "P")
                {
                        print OUT $k."\t".$hash->{$k}->{'Target'}."\n";
                }
                elsif ($hash->{$k}->{'Ass'} eq "S")
                {
                        print OUT $k."\t".$hash_taxa->{$hash->{$k}->{'Target'}}->{'PHYLUM'}."\n";
                }  
                else 
                {
                        print OUT2 $k."\t".$hash->{$k}->{'Target'}."\n";
          	} 
	}
	else 
	{
		print OUT2 $k."\t".$hash->{$k}->{'Target'}."\t".$hash->{$k}->{'identity'}."\n";
	}
}





