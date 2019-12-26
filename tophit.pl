#!/usr/bin/perl

my $file = $ARGV[0]; # Blast Output file to find out the functions of respective hits
open (BLASTM8,"$file") || die "Error: opening\n";
open (OUT,">$file.TopHit") || die "Error: opening TOPHIT\n";

my $blast={};
while(chomp(my $line1=<BLASTM8>)) 
{
	my @arr1=split("\t",$line1);
	my $query=$arr1[0];
	if(!exists($blast->{$query}))
	{
		$blast->{$query}++;
		print OUT "$line1\n"; 
	}
}
close BLASTM8; close OUT1;
