#!/usr/bin/perl
# Author: Ashok K. Sharma

use strict;
use Data::Dumper;

chomp (my $inp=$ARGV[0]);

open(IF,"$inp") or die "Cannot open input file\n"; 
open(OF,">$inp.score") or die "Cannot open input file\n";

##########################################################
my $F={}; my $R={};
####################################################

my $blastoutF = {};my $score;

while(chomp(my $line = <IF>)) 
{
	if ($line =~ /^#/)
	{
		next
	}
        my @arr=split(/\t/,$line);
        #print "start $arr[5]\n";
        if(exists($blastoutF->{$arr[0]}))
        {
                # print "A1 $arr[5]\n";<stdin>;
                if($arr[11] == $score)
                {
                         push (@{$blastoutF->{$arr[0]}}, $line);
                }
		elsif ($arr[11] > $score)
		{
			@{$blastoutF->{$arr[0]}} = ();
			push (@{$blastoutF->{$arr[0]}}, $line);
		}
                next if(!eof);
        }
	else {
        push (@{$blastoutF->{$arr[0]}}, $line);;
        $score=$arr[11];
		}
        #print "SCORE=$score\t$line\n";
}
#print Dumper $blastoutF;

close(IF);

print "blastF array done....printing in progress\n";

foreach my $query (keys %{$blastoutF})
{
        # print "Processing started";
        foreach my $hit (@{$blastoutF->{$query}})
        {
                        print OF $hit."\n";
        
        }
}
  
