#!/usr/bin/perl 
## Author: Ashok K. Sharma
##==>CAZyDB-ec-info.txt.07-20-2017 <== 
#genbank	#family	ec	name
#ABC27701.1	AA10	1.-.-.-	lytic cellulose monooxygenase (HcAA10-2;HCH_00807)
#NP_465990.1	AA10	1.-.-.-	lytic chitin monooxygenase / GlcNAc binding protein A (GbpA;LmLPMO10;Lmo2467)

use strict; 
use Data::Dumper; 

my $infile1 = "/home/gomeza/sharm646/Databases/CaZyDB_130918/CAZyDB-ec-info.txt.07-20-2017";
open(IN1, $infile1) || die $!;

my $cazyhits = {};
while(chomp(my $line = <IN1>))
{
        my @arr = split(/\t/, $line);
        my $id = $arr[0];
        my $EC = $arr[2];
	my $fam = $arr[1];
	$cazyhits->{$id}="$fam\t$EC";
}
close(IN1);
#print Dumper $cazyhits;

my $infile2 = "/home/gomeza/sharm646/Metagenomic/CaZy_Analysis/Total_against_Cazy.BlastOut.besthit";
open (IN2, $infile2) || die $!;
open (OUT, ">$infile2.CazyAnot") || die $!;
while(chomp(my $line1 = <IN2>))
{
#Gene_ID_HU100_S17.genes_35	APK85528.1|GH32	100.000	477	0	0	1	477	1	477	0.0	994	99
#Gene_ID_HU100_S17.genes_124	APK80592.1|GT83	100.000	550	0	0	1	550	1	550	0.0	1116	99
	my @arr1 = split("\t", $line1);
	my @arr2 = split(/\|/, $arr1[1]);
	if (exists ($cazyhits->{$arr2[0]}))
	{
		print OUT $arr1[0]."\t"."$arr2[0]\t".$cazyhits->{$arr2[0]}."\n";
	}
	else
	{
		print OUT $arr1[0]."\t"."$arr2[0]\t"."$arr2[1]\tNA\n";
	}
}

	
