# Author: Ashok K. Sharma

use strict;
use Data::Dumper;

my $shan= {};

my $simp = {};

my $chao = {};
my $Gene = {};

for (my $i= 0; $i <=9 ; $i++)
{
	open (IF, "alpha_rarefaction_4400000_$i.txt") || die "Cannot open $i file\n";
	
	chomp (my $line = <IF>);
	
	while (chomp ($line = <IF>))
	{
		my @arr = split ("\t", $line);
		
		my $sample = shift @arr;
		
		$shan->{$sample}+=$arr[0];
		$simp->{$sample}+=$arr[1];
		$chao->{$sample}+=$arr[2];
		$Gene->{$sample}+=$arr[3];

	}
	
	close IF;
}
print Dumper $shan;
print Dumper $simp;
print Dumper $chao;
print Dumper $Gene;

open (OF, ">Alpha_diversity_mean") || die "Error Out\n";

print OF "\tShannon\tSimpson\tChao1\tObserved Genes\n";

foreach my $val (keys %{$shan})
{
	my $sn = ($shan->{$val})/10;
	
	my $sm = ($simp->{$val})/10;
	
	my $ch = ($chao->{$val})/10;

	my $ge = ($Gene->{$val})/10;
	
	print OF "$val\t$sn\t$sm\t$ch\t$ge\n";
}	

close OF;
