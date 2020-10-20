#!/usr/bin/perl
# Author: Ashok K. Sharma

#==> KEGG_24112010_ko.list <==
#path:ko00010    ko:K00001       ko:E1.1.1.1 ec:1.1.1.1

#==> KEGG_24112010_ko <==
#ENTRY       K00001                      KO

#==> KEGG_24112010_map_title.tab <==
#00010   Glycolysis / Gluconeogenesis

#==> MOUSE_SFB.fa.SFBMouse_00019.pep_pathwaymapping.besthit <==
#MOUSE_SFB|1370722_1371561|-|Exact_GlMet|1620043|Within|complete cno:NT01CX_1146 45.49   277     148     1       2       275     3       279     2e-67    258


use strict;
use Data::Dumper;

#print "enter the best hit file\n";
chomp (my $infile4=<$ARGV[0]>);

#my $path = "/home2/darshan/KEGG_Analysis_scripts";
my $infile1 = "/home/gomeza/sharm646/Databases/KEGG_Database/KEGG_2014_ko.list";
open(IN1, $infile1) || die $!;

my $pathway = {};
while(chomp(my $line = <IN1>))
{
	my @arr = split(/\t/, $line);
	$arr[0]=~/^\w+:(ko\d+)/;
	my $path_id = $1;
	$arr[1]=~/^ko:(K\d+)/;
	my $ko_id = $1;
	push(@{$pathway->{$ko_id}} , $path_id);
}
close(IN1);

my $infile2 = "/home/gomeza/sharm646/Databases/KEGG_Database/KEGG_2014_ko";
open(IN2, $infile2) || die $!;

my $ko = {};
my $desc;
my $id;

while(chomp(my $line = <IN2>))
{
	if($line=~/^ENTRY\s+(K\d+)/)
	{
		$id = $1;
	}
	
	if($line=~/^DEFINITION\s+(.*)/)
	{
		$desc = $1;
	}

	if($line=~/^GENES\s+(.*)$/)
	{
		my @rows;
		push(@rows, $1);	
		while(chomp(my $line=<IN2>))	
		{
			if($line=~/^(\/\/\/|\w+)/ || eof)
			{
				last;
			}
			push(@rows,$line);
		}
		my $prefix;
		my @genes;
		foreach my $row (@rows)
		{
			$row=~s/\(.*?\)//g;
			$row=~s/^\s+|\s+$//;
			my @arr=split(/\s+/, $row);
			if($arr[0]=~/:$/)
			{
				$prefix = lc(shift @arr);
			}
			
			push(@genes, map {$_=$prefix.$_} @arr);
		} 
		
		foreach my $gene (@genes)
		{
			my $na = ["NA"];
			my $pathway_id = exists($pathway->{$id}) ?  $pathway->{$id} : $na;
			$ko->{$gene}=[$id, $desc, $pathway_id];
		}
	}
}
close(IN2);

my $infile5 = "/home/gomeza/sharm646/Databases/KEGG_Database/KEGG_GENOPEP_200314.single";
open (IN5, $infile5) || die $!;

my $hash_length = {};
my $annot;
my $seq_length;
while (chomp (my $line = <IN5>))
{
	if ($line =~ /^>/)
	{
		($annot) = $line =~ />(\S+:\S+)/; 
	}
	else {
		my @seq = split ('', $line);
		$seq_length = scalar (@seq);
	     }
	$hash_length->{$annot}=$seq_length;
}

my $infile3 = "/home/gomeza/sharm646/Databases/KEGG_Database/KEGG_map_title.tab_2014";
open(IN3, $infile3) || die $!;

my $title={};
while(chomp(my $line = <IN3>))
{
	my @arr = split(/\t/, $line);
	$title->{"ko".$arr[0]} = $arr[1];
}
close(IN3);

#my $infile4 = "RAT_SFB_Start_Corrected_SEQ.pep_pathwaymapping.besthit";
#my $infile4 = "RAT_SFB_FinalORFs_AfterManCheck.faa_pathwaymapping.bestID_on_Coverage";
my $infile4 = <$ARGV[0]>;
open(IN4, $infile4) || die $!;

my $group;

my $HASH = {};
open(OUT1, ">$infile4.KO") || die $!;
open (REJECT, ">$infile4.REJECT") || die $!;
while(chomp(my $line = <IN4>))
{
	my @arr = split(/\t/, $line);
	if (($arr[-2] >= 60) && ($arr[-3] <= 0.000001))
	{
	if (!exists($HASH->{$arr[0]}))
	{
	if(exists($ko->{$arr[1]}))
	{
		my $kos = join(';', @{$ko->{$arr[1]}->[2]});
		my $titles = join(';', map {$title->{$_}}  @{$ko->{$arr[1]}->[2]});
		foreach my $k ( @{$ko->{$arr[1]}->[2]})
		{
			next if($k eq "NA");
			if(!grep(/^$ko->{$arr[1]}->[0]$/,@{$group->{$k}}))
			{
				push(@{$group->{$k}}, $ko->{$arr[1]}->[0]);
			}
		}
		my $ratio = $arr[3]/$hash_length->{$arr[1]};
		print OUT1 $arr[0] ."\t". $arr[1] ."\t". $ko->{$arr[1]}->[0] ."\t". $ko->{$arr[1]}->[1] ."\t". $kos . "\t". (($titles) ? $titles: "NA") ."\t".$ratio."\n";
		$HASH->{$arr[0]}++;
	}
	else
	{
		print OUT1 $arr[0] ."\t" .$arr[1] ."\t". "NA"  ."\t". "NA" ."\t". "NA" ."\t". "NA"  ."\n";
		$HASH->{$arr[0]}++;
	}
	}
	else {
		next
		}
	}
	else 
	{
		print REJECT $line."\n";
	}
}
close(IN4);
close(OUT1);

open(OUT2, ">$infile4.map") || die $!;
foreach my $e (keys %{$group})
{
	print OUT2 $e ."+". join('+', map{ $_="ko:$_"} @{$group->{$e}}) ."\n";
}
close(OUT2);
exit;
