#!/usr/bin/perl -w
use Data::Dumper;
use strict;
my $input = <$ARGV[0]>;
open (IN, $input) || die "Cannot open the file\n"; # Input a fna Total protein file file
my $hash={};

while (chomp(my $line = <IN>))
{
        if ($line =~ />/)
        {
		my $ID = $line;
		chomp(my $seq = <IN>);
		$hash->{$ID}=$seq;
	}
}
#print Dumper $hash;

my $input1 = <$ARGV[1]>;
open (IN1, $input1) || die "Cannot opent he file\n"; # Input a total fna clustered file
open (OUT, ">Combined_Total_Proteins.NRset.faa") || die "Cannot open the file\n";
while (chomp(my $line1 = <IN1>))
{
	if ($line1 =~ />/)
	{
                my $ID1 = $line1;
                $ID1 =~ s/genes_/proteins_/g;
		#print $ID1; <stdin>;
		if (exists($hash->{$ID1}))
                {
			my $ID2 = $ID1;
			$ID2 =~ s/proteins_/genes_/g;
                        print OUT $ID2."\n".$hash->{$ID1}."\n";
                }
	}
}

