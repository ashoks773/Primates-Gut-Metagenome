#!/usr/bin/perl
# Author: Ashok K. Sharma

use strict;

my $infile=$ARGV[0];
open (BLASTM8,$infile) || die "Error: opening $infile\n";
open (OUT,">$infile.besthit") || die "Error: opening $infile.besthit\n";
while(chomp(my $line=<BLASTM8>)) 
{
	my @arr=split(/\t/,$line);
	my $iden=$arr[2];
	my $qcov=$arr[12];
	my $score=$arr[11];
	my $eval=$arr[10];
	my @arr1 = split (/\-/, $eval);
	#print $arr1[1]; <stdin>;	
	if($qcov >= 80 && $arr1[1] >= 6 && $iden >= 50 && $score >= 60)
	{
	#	print "$line\n"; <stdin>;
		print OUT "$line\n"; 
	}
	if($qcov >= 80 && $eval==0 && $iden >= 50 && $score >= 60)
	{
	#	print "$line\n"; <stdin>;
                print OUT "$line\n";
	}
}
close(BLASTM8);
exit;

