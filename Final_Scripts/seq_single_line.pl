#!usr/bin/perl
##This program reads a GenBank file and converts it into a fasta format
#file with sequenc in a single line

print "Enter the name of input file\n";
chomp($file=<stdin>);
$seq="";
open (IF,"$file") || die "cannot open the input file\n";
open (OF,">$file.single") || die "cannot open the $! file\n";

while(chomp($line=<IF>))
{
  if($line =~/^>/)
   {
     print OF  "$seq\n";
     $seq="";$annot="";
     $annot=$line;
     print OF "$annot\n";
   }
 else
   {
     $seq.=$line;
   }
}
print OF "$seq\n";
exit;
