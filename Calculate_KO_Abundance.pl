#!/usr/bin/perl -w
#
#This script Calculates the KO and EC Abundance of Genes from the Gene Table

use strict;

my $infile = <$ARGV[0]>;

open (IN, $infile) || die;

my $hash = {};
my $hash2 = {};
my $hash3={};
while (chomp(my $line = <IN>))
{
	my @x = split ("\t", $line);
	if ($x[1] ne "NA")
	{
	my $KO = $x[2];
#	my $annot = $x[3];
#	$hash3->{$KO}->{'KO'}=$annot;
	my $gene_ID = $x[0];
	$hash->{$gene_ID}->{'KO'}=$KO;
#	$hash->{$KO}->{'annot'}=$annot;
	}
#        (my $annot) = $x[3] =~ /(.*)\[/;
        my $EC = $x[7];
	my $gene_ID = $x[0];
#	$hash3->{$EC}->{'EC'}=$annot;
#        if ($EC ne "")
#        {
	$hash->{$gene_ID}->{'EC'}=$EC;
#        }
}


my $Gene_abundance_file = <$ARGV[1]>;
open (GENE_ABUNDANCE, $Gene_abundance_file) || die "Cannot open the file\n";

my @v;

my $EC_hash = {};
while (chomp(my $line2 = <GENE_ABUNDANCE>))
{
        if ($line2 =~ /^Genes/)
        {
                @v = split ("\t", $line2);
                shift (@v);
                next
        }
        else {
                my @t = split ("\t", $line2);
                my $GENE_id = $t[0];
                if (exists($hash->{$GENE_id}->{'KO'}))
                {
                        shift(@t);
                        foreach my $v(@v)
                        {
				$hash2->{$hash->{$GENE_id}->{'KO'}}->{$v}+=$t[0];
				if (exists($hash->{$GENE_id}->{'EC'}))
				{
				$EC_hash->{$hash->{$GENE_id}->{'EC'}}->{$v}+=$t[0];
				}
                                shift (@t);
                        }
                }

}
}
open (OUT, ">KO_ABUNDANCE") || die "Cannot open the file\n";

open (OUT2, ">EC_ABUNDANCE") || die "Cannot open the file\n";



unshift (@v, "KO_ID");

my $Header = join ("\t", @v);

print OUT $Header."\n";

print OUT2 $Header."\n";

shift (@v);

foreach my $r  (keys %{$hash2})
{
        print OUT $r;
        foreach my $v (@v)
        {
                print OUT "\t".$hash2->{$r}->{$v};
        }
        print OUT "\n";
}

foreach my $c (sort keys %{$EC_hash})
{
        print OUT2 $c;
        foreach my $v (@v)
        {
                print OUT2 "\t".$EC_hash->{$c}->{$v};
        }
	print OUT2 "\n";
}
