#!/usr/bin/perl 
# Author: Ashok K. Sharma

use strict; 
use warnings; 
use Data::Dumper;

my $GENES = <$ARGV[0]>; #Total NR Gene set
open (IN, $GENES) || die "Cannot open the gene file\n";

my $Length_hash = {};
while (chomp(my $line = <IN>))
{
        if ($line =~ />/)
        {
                (my $I) = $line =~ />(\S*)/;
                chomp ($line = <IN>);
                my @SEQ = split ("", $line);                 ## Get the gene lengths into a hash
                my $SEQ_LENGTH = scalar(@SEQ);
                $Length_hash->{$I}=$SEQ_LENGTH;
        }
}

#print Dumper $Length_hash;
my $genecounts = <$ARGV[1]>; #Gene counts in all the samples
open (IN1, $genecounts) || die "Cannot open the gene countsfile\n"; 
open (OUT, ">GENE_ABUNDANCE_TABLE");
while (chomp(my $line1 = <IN1>)) 
{
	my @arr = split ("\t", $line1);
	my $gname = shift (@arr);
	if (exists $Length_hash->{$gname})
	{
		#print "$Length_hash->{$gname}"; <stdin>;
		print OUT "$gname";
		foreach my $k (@arr)
		{
		#	print "$k\n"; <stdin>; 	print "$Length_hash->{$gname}"; <stdin>;
			my $Fraction = $k/$Length_hash->{$gname}; # print "\n$Fraction\n"; <stdin>;
			my $rounded = sprintf ("%.4f", $Fraction);
			print OUT "\t$rounded";
		}
		print OUT "\n";
	}
}

		
