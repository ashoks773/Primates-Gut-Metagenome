#!/usr/bin/perl -w

use strict;

my $input = <$ARGV[0]>;
open (IN, $input) || die "Cannot open the file\n";

open (OUT, ">$input.filtered") || die "Cannot open the file\n";

chomp(my $line = <IN>);
print OUT $line."\n";

while (chomp(my $line = <IN>))
{
	my $N = 0;
	my @x = split ("\t", $line);
	my $ID = shift(@x);
	foreach my $x (@x)
	{
		if ($x > 0.001)
		{
			$N++
		}
	}
	if ($N >= 3)
	{
		print OUT $line."\n";
	}
}

