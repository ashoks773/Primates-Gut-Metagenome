#!/usr/bin/perl -w
## Author: Ashok K. Sharma
use strict;
use Data::Dumper;

my $hash= {};

my $file=$ARGV[0];
open (IF, $file) || die "Can't open a input file \n";
while (chomp (my $line = <IF>))
{
#Samples Group   Group1  Group2
#2009_106        Gorilla Dry     Habituated

	my @arr = split ("\t", $line);
	$hash->{$arr[0]}=$line;
}
#print Dumper $hash;

my $file1=$ARGV[1];
open (IF1, $file1) || die "Can't open the second file \n";
open (OF, ">$file1.out") || die "Can't open the output file \n";
while (chomp (my $line1 = <IF1>))
{
	my @arr1 = split ("\t", $line1);
	if (exists ($hash->{$arr1[0]}))
	{
		print OF "$line1\t$hash->{$arr1[0]}\n";
	}
}






