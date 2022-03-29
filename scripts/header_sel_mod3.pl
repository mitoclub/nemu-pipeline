#!/usr/bin/perl -w
$inxx=0; $outxx=0;
$in=0; $out=0;
%ingroup=();
%outgroup=();
%ingroupnames=();
%outgroupnames=();
@namelist=();
$name=$ARGV[1];
$name=~s/_/ /g;
open (I, "$ARGV[0]");
while (<I>)
{
#fasta header
if ($_=~/^>(.*?)\n/)
    {
    $tmp=$1;
    if ($tmp=~/$name/)
	{
	$inxx++;
	$replaced_name="RN_$inxx";
	$ingroupnames{$replaced_name}=$tmp;
	$in=1; $out=0;
	}
    else
	{
	$outxx++;
	$replaced_name="OUTGRP_$outxx";
	$outgroupnames{$replaced_name}=$tmp;
	$in=0; $out=1;
	}
    }
#fasta sequence
else
    {
    push @namelist, $replaced_name;
    if ($in==1 && $out==0) {$ingroup{$replaced_name}.=$_}
    if ($in==0 && $out==1) {$outgroup{$replaced_name}.=$_}
    }
}
close I;

if ($namelist[$#namelist]=~/OUTGRP_/)
{
$outgrp="";
for ($j=$#namelist-1;$j>-1;$j--)
    {
    if ($namelist[$j]=~/RN_/ and $namelist[$j+1]=~/OUTGRP_/) {$outgrp=$namelist[$j+1]; last}
    }
}
elsif ($outxx>0)
{
for ($j=$#namelist;$j>-1;$j--)
    {
    if ($namelist[$j]=~/OUTGRP_/) {$outgrp=$namelist[$j]; last}
    }
}
else
{
$outgrp=$namelist[$#namelist];
$outgroup{$outgrp}=$ingroup{$outgrp};
$outgroupnames{$outgrp}=$ingroupnames{$outgrp};
}

print ">OUTGRP\n$outgroup{$outgrp}";
print STDERR "OUTGRP\n$outgroupnames{$outgrp}\n";
foreach $in (keys %ingroupnames)
{
print ">$in\n$ingroup{$in}";
print STDERR "$in\n$ingroupnames{$in}\n";
}
