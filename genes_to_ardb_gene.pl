#!/usr/bin/perl 
#
##==> ARDB Genes <== 
#YP_637442	Aac2I	164756	NULL	C_NC_008146.1|305538|306161
#ABG06386	Aac2I	164756	NULL	C_CP000384.1|305538|306161

use strict; 
use Data::Dumper; 

my $infile = "/home/gomeza/sharm646/Databases/ardbAnno1.0/tabs/abrg.tab"; 
open(IN1, $infile) || die $!;
my $ardbhits = {};
while(chomp(my $line = <IN1>))
{
        my @arr = split(/\t/, $line);
        my $id = $arr[0];
        my $ardb_id = $arr[1];
	#my $fam = $arr[1];
	$ardbhits->{$id}="$ardb_id";
}
close(IN1);
#print Dumper $ardbhits;

my $infile2 = "/home/gomeza/sharm646/Metagenomic/Antibiotic_Resistance_ARDB/Total_against_ARDB.BlastOut";
open (IN2, $infile2) || die $!;
open (OUT, ">$infile2.ARDB_anot") || die $!;
while(chomp(my $line1 = <IN2>))
{
#Gene_ID_HU100_S17.genes_126	sp|B6I7J8|	100.000	660	0	0	1	660	1	660	0.0	1368	99	100.00
#Gene_ID_HU100_S17.genes_186	dbj|BAG77975|	100.000	396	0	0	1	396	1	396	0.0	788	99	100.00
	my @arr1 = split("\t", $line1);
	my @arr2 = split(/\|/, $arr1[1]);
	if (exists ($ardbhits->{$arr2[1]}))
	{
		print OUT $arr1[0]."\t"."$arr2[1]\t".$ardbhits->{$arr2[1]}."\n";
	}

	#else
	#{
	#	print OUT $arr1[0]."\t"."$arr2[0]\t"."$arr2[1]\tNA\n";
	#}
}

	
