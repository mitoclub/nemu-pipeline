#!/usr/bin/perl -w

%ML=();
open (I, "$ARGV[0]");
while(<I>)
{
chomp;
@w=split(/\s+/);
$nm=shift @w;
@I=split(//, shift @w); 
@{$ML{$nm}}=@I
}
close I;

%MP=();
open (I, "$ARGV[1]");
while(<I>)
{
chomp;
@w=split(/\s+/);
$nm=shift @w;
@I=split(//, shift @w);
@{$MP{$nm}}=@I
}
close I;

%deep=();
$ln=0;
open (I, "$ARGV[2]");
while(<I>)
{
$ln++;
chomp;
if ($ln==1)
    {
    @nms=split(/\s+/, $_);
    shift @nms;
    }
elsif ($_=~/^OUTGRP/)
    {
    @w=split(/\s+/, $_);
    shift @w;
    for ($j=0; $j<$#w+1; $j++)
	{
	$deep{$nms[$j]}=$w[$j];
	}
    }
}
close I;

open (O, ">$ARGV[3]");
print "Tree node\tDeep\tDissimilarity\tPositions\n";
foreach $nm (sort {$a cmp $b} keys %ML)
{
@consensus=();
@positions=();
$diss=0;
for ($i=0; $i<$#{$ML{$nm}}+1; $i++)
    {
    if ($ML{$nm}[$i] eq $MP{$nm}[$i]) {push @consensus, $ML{$nm}[$i]}
    else {push @consensus, "N"; $diss++; push @positions, $i+1}
    }
$positions=join(";", @positions);
if ($diss>0) {print "$nm\t$deep{$nm}\t$diss\t$positions\n"}
$consensus=join("", @consensus);
print O "$nm $consensus\n";
}
close O;
