#!/usr/bin/perl -w
# Author: Ashok K. Sharma

use strict;

my $input = <$ARGV[0]>;
open (IN, $input) || die "Cannot open the file\n";

open (OUT, ">$input.Filtered") || die;


my $ID;

while (chomp(my $line = <IN>))
{
LINE:	if ($line =~ />/)
	{
		$ID = $line;
#		print $line;<stdin>;
	}
	else {
		my $seq;
		until ($line =~ />/)
		{
			$seq .= $line;
#			print $line;<stdin>;
			last if eof(IN);
			chomp($line = <IN>);
		}
#		print $seq;<stdin>;
		my @S = split ("", $seq);
		my $N = $seq =~ s/N//g;
		if ($N == 0)
		{
		if (scalar(@S) >= 1000)
		{
			print OUT $ID."\n".$seq."\n";
		}
		}
		last if eof(IN);
		goto LINE;
		}
}

