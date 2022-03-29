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
if ($seqNumber>=4)
{
open (Y, ">yes.O");
print Y "Yes $seqNumber!\n";
close Y;
}

print "$outgrp";
foreach $seq (keys %alig)
{
print ">$alig{$seq}\n$seq\n";
}
