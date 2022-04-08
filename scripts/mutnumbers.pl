#!/usr/bin/perl -w

%terminals=();
$seq="";
open (I, "$ARGV[0]");
while(<I>)
{
chomp;
if ($_=~/>(\w+)$/) 
    {
    $nm=$1;
    if ($seq=~/\w/) {$seq=~s/\s+//g; @seq=split(//, $seq); @{$terminals{$name}}=@seq; $seq="";}
    $name=$nm;
    }
else 
    {
    $seq.=$_;
    }
}
if ($seq=~/\w/) {$seq=~s/\s+//g; @seq=split(//, $seq); @{$terminals{$name}}=@seq; $seq=""}
$max=$#seq+1;
close I;

delete $terminals{OUTGRP};

$minMut=$maxMut=0;
print "Position\t# of variants\tNucleotides\n";
for ($i=0;$i<$max;$i++)
{
%variants=();
foreach $name (sort {$a cmp $b} keys %terminals)
    {
    if ($terminals{$name}[$i]=~/[ATGC]/) {$variants{$terminals{$name}[$i]}+=1}
    }
$variantsN=keys %variants;
$s=0;
$variantsS="";
if ($variantsN>1)
{
foreach $v (sort {$variants{$a}<=>$variants{$b}} keys %variants)
    {
    $s++;
    $mass+=$variants{$v};
    unless ($s==$variantsN) {$maxMut+=$variants{$v}}
    $variantsS.=";$v($variants{$v})";
    }
$variantsS=~s/^;//;
print "$i\t$variantsN\t$variantsS\n"; 
$minMut++;
}
}

if ($max > 0) {
    $p=$minMut/$max;
    print "Variable positions\t$minMut\t$max\t$p\n";
}
if ($mass > 0) {
    $p=$maxMut/$mass;
    print "Maximum number of substitutions\t$maxMut\t$mass\t$p\n";
}
