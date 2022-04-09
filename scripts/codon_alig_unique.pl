#!/usr/bin/perl -w

%alig=();
$seq="";
$outgrp="";
open (I, "$ARGV[0]");
while(<I>)
{
chomp;
if ($_=~/^>/) 
    {
    $nm=$_; 
    $nm=~s/>//;
    if ($seq=~/\w/) 
	{
	if ($name=~/RN_\d+/) {$seq=~s/\s+//g; $alig{$seq}=$name}
	elsif ($name=~/OUTGRP/) {$seq=~s/\s+//g; $outgrp=">$name\n$seq\n"}
	$seq="";
	}
    $name=$nm;
    }
else
    {$seq.=$_}
}
close I;
if ($seq=~/\w/)
    {
    if ($name=~/RN_\d+/) {$seq=~s/\s+//g; $alig{$seq}=$name}
    elsif ($name=~/OUTGRP/) {$seq=~s/\s+//g; $outgrp=">$name\n$seq\n"}
    }

$seqNumber=keys %alig;
if ($seqNumber>=4) {
    open (Y, ">report_yes.txt");
    print Y "Seq number is $seqNumber!\n";
    close Y;
} else {
    open (N, ">report_no.txt");
    print N "Seq number is too low - $seqNumber!\nBreak execution!\n";
    close N;
}

print "$outgrp";
foreach $seq (keys %alig)
{
print ">$alig{$seq}\n$seq\n";
}
