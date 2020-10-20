#!/usr/bin/perl -w
# Author: Ashok K. Sharma

use strict;
use Data::Dumper;
my $input = <$ARGV[0]>;
open (IN, $input) || die "Canno topent he file\n";
#CP002006.faa	Prevotella ruminicola 23
#CP002122.faa	Prevotella melaninogenica ATCC 25845

my $hash={};
while (chomp(my $line = <IN>))
{
	my @x = split ("\t", $line);
	$hash->{$x[0]}=$x[1];
}
#print Dumper $hash;

my $input2 = <$ARGV[1]>;
open (IN2, $input2) || die "can't open the file\n";
#>AFPX00000000.faa
#MRNAVRKFRN ISDKEEAVKL YPTLQKMLDK LAKVNIIHKN KAANLKSGLC HHIATLGSHT
open (OUT, ">$input2.Names.aln") || die "Cannot open the file\n";
while (chomp(my $line2 = <IN2>))
{
	if ($line2 =~ /^>/)
	{
		$line2 =~ s/>//g;
		if (exists($hash->{$line2}))
		{
			print OUT ">$hash->{$line2}\n";
		}
		if (!exists($hash->{$line2}))
		{
			print OUT ">$line2\n";
		}
	}
	else {
		print OUT "$line2\n";
	}
}
