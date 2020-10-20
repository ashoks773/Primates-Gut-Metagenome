#!/usr/bin/perl -w
# Author: Ashok K. Sharma

use strict;

#my $input = <$ARGV[0]>;
my $input = "/home/gomeza/sharm646/Databases/KEGG_Database/KEGG_2014_ko.list";
open (IN, $input) || die "Cannot open the file\n";

#path:ko00010	ko:K00001	ko:E1.1.1.1 ec:1.1.1.1
#path:ko00010	ko:K00002	ko:AKR1A1 ec:1.1.1.2

my $hash={};
my $hash_annot = {};
my $hash_pathcount={};

while (chomp(my $line = <IN>))
{
	my @x = split ("\t", $line);
	(my $KO_ID) = $x[1] =~ /ko:(\S+)/;
	(my $PATH_ID) = $x[0] =~ /path:ko(\S+)/;
	if (defined($KO_ID))
	{
		push (@{$hash->{$KO_ID}}, $PATH_ID);
	}
}

#my $input2 = <$ARGV[1]>;
my $input2 = "/home/gomeza/sharm646/Databases/KEGG_Database/KEGG_map_title.tab_2014";
open (IN2, $input2) || die "Cannot open the file\n";

#00010	Glycolysis / Gluconeogenesis
#00020	Citrate cycle (TCA cycle)

while (chomp(my $line = <IN2>))
{
	my @y = split ("\t", $line);
	my $path = $y[0];
	my $annot = $y[1];
	$hash_annot->{$path}=$annot;
}

open (OUT, ">PATHWAY_ABUNDANCE_TABLE") || die "Cannot open the file\n";

my $input3 = <$ARGV[0]>; # KO Abundacne file
open (IN3, $input3) || die "Cannot open the file\n";

#EC_ID	AB_0	AB_C	AB_M	AB_W	T_0	T_C	T_M	T_W	U_0	U_C	U_M	U_W	V_0	V_C	V_M	V_W	Y_0	Y_C	Y_M	Y_W	Z_0	Z_C	Z_M	Z_W
#K13002	0.001506	0.025387

open (OUT2, ">CHECK_CALCULATION") || die "Cannot open the file\n";

my @v = ();
while (chomp(my $line = <IN3>))
{
	if ($line =~ /^KO_ID/)
	{
		@v = split ("\t", $line);
		shift (@v);
	}	
	else {
		my @x = split ("\t", $line);
		my $KO_ID = shift (@x);
		if (exists($hash->{$KO_ID}))
		{
			my $n = scalar (@{$hash->{$KO_ID}});
			foreach my $h (@{$hash->{$KO_ID}})
			{
				print OUT2 $KO_ID."\t".$h."\t".$x[0]."\t".$n."\n";
			}
			foreach my $v (@v)
			{
				foreach my $p (@{$hash->{$KO_ID}})
				{
					$hash_pathcount->{$p}->{$v}+=$x[0]/$n;
				}
				shift (@x);
			}
		}
	}
}

my $header = join ("\t", @v);
print OUT "Pathway"."\t".$header."\n";

foreach my $k (keys %{$hash_pathcount})
{
	print OUT $hash_annot->{$k};
	foreach my $v (@v)
	{
		print OUT "\t".$hash_pathcount->{$k}->{$v};
	}
	print OUT "\n";
	
}
