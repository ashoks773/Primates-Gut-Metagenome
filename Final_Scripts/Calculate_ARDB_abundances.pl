#!/usr/bin/perl -w
#This script Calculates the ARDB Gene abundance from the Gene Table
# Author: Ashok K. Sharma

use strict;
use Data::Dumper;

my $infile = <$ARGV[0]>;
open (IN, $infile) || die;
my $hash = {};
my $hash2 = {};
my $hash3={};
while (chomp(my $line = <IN>))
{
	my @x = split ("\t", $line);
	my $ardbgene = $x[2];
	my $gene_ID = $x[0];
	$hash->{$gene_ID}->{'ARDBGENE'}=$ardbgene;
        #my $EC = $x[3];
	#my $gene_ID = $x[0];
	#$hash->{$gene_ID}->{'EC'}=$EC;
	#my $genebank = $x[1];
	#$hash->{$gene_ID}->{'GeneBank'}=$genebank;
}
#print Dumper $hash;

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
                if (exists($hash->{$GENE_id}->{'ARDBGENE'}))
                {
                        shift(@t);
                        foreach my $v(@v)
                        {
				$hash2->{$hash->{$GENE_id}->{'ARDBGENE'}}->{$v}+=$t[0];
                                shift (@t);
                        }
                }

}
}
open (OUT, ">ARDB_GENE_ABUNDANCE") || die "Cannot open the file\n";
unshift (@v, "ARDB_Gene");
my $Header = join ("\t", @v);

print OUT $Header."\n";
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

