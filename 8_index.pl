#!/usr/bin/perl -w
use strict;

my $N = 0;
foreach my $arg(@ARGV)
{
	$N++;
	(my $name) = $arg =~ /^(\S+)\.fna/;
	open (IN, $arg) || die "Cannot open the file\n";
	my $n = 1;
	open (OUT, ">".$name."_indexed_GENES.fna") || die "Cannot open the output file\n";
	while (chomp (my $line = <IN>))
	{
		if ($line =~ />/)
		{
			print OUT ">Gene_ID_".$name."_".$n."\n";
			$n++;
		}
		else {
			print OUT $line."\n";
			}
	}
}
close (IN);
close (OUT);
