#!/usr/bin/perl -w
# USAGE: codon_alig_unique.pl seqs.fasta  1>seqs_unique.fasta

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

print "$outgrp";
foreach $seq (keys %alig)
{
print ">$alig{$seq}\n$seq\n";
}
