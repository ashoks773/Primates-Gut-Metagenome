#!/usr/bin/perl -w
#
use strict;

#This script takes as input the Genome ID file, the Genes file and the Blast output file and gives annotation to the Blast hits based on identity and aligned length

#Achromobacter_xylosoxidans	gi|311103224|ref|NC_014640.1|	gi|311109684|ref|NC_014641.1|	gi|311109789|ref|NC_014642.1|
#Achromobacter_xylosoxidans	gi|568125122|ref|NC_023061.1|


my $Genome_ID = <$ARGV[0]>;
open (GEN_ID, $Genome_ID) || die "Cannot open the file\n";

my $hash={};
my $hash_length = {};
my $hash_hit = {};

while (chomp (my $Ref_line = <GEN_ID>))
{
	my @x = split ("\t", $Ref_line);
	my $Genome = shift(@x);
	foreach my $x (@x)
	{
		$hash->{$x}=$Genome;
	}
}



#>gene_id_18_156096
#AGGGAGACAGCGGGCAAGCTGCTCGATATGGGGTTCTATATCTCGTTTTCCGGTATCGTCACTTTTCGTAACGCCGAGCAGCTGCGTGATGCGGCGCGCTACGTGCCGCTGGATCGTCTGCTGGTTGAAACCGACTCGCCGTACCTGGCGCCGGTCCCCCACCGCGGCAAAGAGAACCAGCCGGCGCTGGTCCGTGACGTTGCCGAG
#TATATGGCGGTACTGAAGGGGGTGTCGCTTGAGCAGCTGGCTCAGCAGACCACCGACAATTTTGCCACCCTGTTCCACATCGACCCCGCCCGTTTGCAACCTGCTTGA
#>gene_id_15_137813
#AACGACTTTGAGGTCGGCCTGCGCCATGACCTGCCCATCGTCCGCGTCTTCACCTACGACGGTCGCATGACCGGCGCCGCCGACAAGGCTGCTGCCGACGCTCTGTTTGCTGCAGGCAAGAACACCGTCAACGAGCCTCATGTGCTG

my $Genes_file = <$ARGV[1]>;
open (GENE_FILE, $Genes_file) || die "Cannot open the file\n";

while (chomp (my $Genes_line = <GENE_FILE>))
{
	if ($Genes_line =~ />/)
	{
		(my $gene_id) = $Genes_line =~ />(\S+)/;
		chomp ($Genes_line = <GENE_FILE>);
		my @l = split ("", $Genes_line);
		my $length = scalar(@l);
		$hash_length->{$gene_id}=$length;
	}
}


#411489.CLOL250_01016   DS480321        100.00  60      0       0       1       60      71227   71286   7e-25    119
###411489.CLOL250_01016   AAYW02000007    100.00  60      0       0       1       60      71227   71286   7e-25    119

my $hits_file = <$ARGV[2]>;
open (HITS, $hits_file) || die "Cannot open the file\n";

open (OUT, ">Hits_Annotation.out") || die "Cannot open the file\n";
open (OUT3, ">Rejected_E_Value.out") || die "Cannot open the file\n";

while (chomp (my $hits_line = <HITS>))
{
	my @hit = split ("\t", $hits_line);
	if (scalar(@hit) > 12)
	{
		next
	}
	else {
	my $query = $hit[0];
	my $target = $hit[1];
	my $identity = $hit[2];
	my $aligned_length = $hit[3];
	my $e_value = $hit[-2];
	my $score = $hit[-1];
	if ($e_value < 0.000001)
	{
		push (@{$hash_hit->{$query}->{'identity'}},$identity);
		push (@{$hash_hit->{$query}->{'AL'}}, $aligned_length);
		push (@{$hash_hit->{$query}->{'TARGET'}}, $target);
	}
	else {
		print OUT3 $hits_line."\n";
		}
}
}



open (OUT1, ">Rejected_hits") || die "Cannot open the file\n";
open (OUT2, ">Unannotated_hits") || die "Cannot open the file\n";

foreach my $k (keys(%{$hash_hit}))
{
	my $max = &max(@{$hash_hit->{$k}->{'identity'}});
	my @removed = splice(@{$hash_hit->{$k}->{'identity'}}, $max);
	splice (@{$hash_hit->{$k}->{'AL'}}, $max);
	splice (@{$hash_hit->{$k}->{'TARGET'}}, $max);
	#my $m = $hash_hit->{$k}->{'identity'}->[0];
	#my $al = $hash_hit->{$k}->{'AL'}->[0];
	if (exists($hash_length->{$k}))
	{
	#my $aligned_percent = ($al/$hash_length->{$k})*100;
	for (my $m = 0; $m < scalar(@{$hash_hit->{$k}->{'TARGET'}}); $m++)
	{
	my $al = $hash_hit->{$k}->{'AL'}->[$m];
	my $identity = $hash_hit->{$k}->{'identity'}->[$m];
	my $aligned_percent = ($al/$hash_length->{$k})*100;
	if ($aligned_percent >= 80)
	{ 
		my $j = $hash_hit->{$k}->{'TARGET'}->[$m];
		if (exists($hash->{$j}))
		{		
			print OUT $k."\t".$hash->{$j}."\t".$identity."\n";
		}
		else {	
			print OUT2 $k."\t".$j."\t".$identity."\n";
			}
	}
	else {	
		print OUT1 $k."\t".$hash_hit->{$k}->{'TARGET'}->[$m]."\t".$hash_hit->{$k}->{'identity'}->[$m]."\t".$hash_hit->{$k}->{'AL'}->[$m]."\t".$aligned_percent."\n";
		}
	}
	}
}



sub max 
{
	my ($max_so_far) = shift (@_);
	my $N=1;
	if (scalar(@_) > 0)
	{
		foreach (@_)
		{
			if ($_ == $max_so_far)
			{
				$N++;
			}
			elsif ($_ > $max_so_far)
			{
				next;
			}
		}
	}
	else 
	{
		$N = 1
	}
	$N;
}


		


 

