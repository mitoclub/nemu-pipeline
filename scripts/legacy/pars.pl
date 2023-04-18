#!/usr/bin/perl -w
@T=();
open (T, "$ARGV[0]");
@T=<T>;
close T;
$W=$T=$T[0];
$T=~s/\)\d+/\)/g;
open (O, ">parsimony.tmp.tre"); print O "$T"; close O;
%branchhash=();
%terminals=();
while ($W=~/\(([^\(\),]+),([^\(\),]+)\)(\d+)/)
{
$branch0=$1;
$branch1=$2;
$ancestor=$3;
if ($branch0=~/[A-Za-z]/) {$terminals{$branch0}=1}
if ($branch1=~/[A-Za-z]/) {$terminals{$branch1}=1}
$branchhash{$branch0}=$ancestor;
$branchhash{$branch1}=$ancestor;
$W=~s/\($branch0,$branch1\)$ancestor/$ancestor/;
#print "$branch0 - $branchhash{$branch0} | $branch1 - $branchhash{$branch1} -> $W\n";
}

foreach $br (keys %branchhash) {print "$br <- $branchhash{$br}\n"}

system ("seqret -osformat2 phylip $ARGV[1] parsimony.tmp.alig");
if (-e "outfile" and -e "outtree") {system ("rm outfile outtree")}
system ("echo -e \"parsimony.tmp.alig\nU\n5\n.\n3\nY\nparsimony.tmp.tre\n\" | ~/bin/phylip/dnapars");


%branchhashP=();
$ln=0;
open (PR, "outfile");
while (<PR>)
{
$ln++;
chomp;
if ($ln>7)
    {
    $_=~s/\s+$//;
    $_=~s/\?|\.|\*/X/g;
    while ($_=~/(\S+) (\S+)/) {$A=$1; $B=$2; $_=~s/$A $B/$A$B/;}
    $str=" $_";
    @str=split(/\s+/, $str);
    if ($#str==4)
	{
	$branchhashP{$str[2]}=$str[1];
	}
    }
}
close PR;
%branchhashPARS=();
$branchhashP=keys %branchhashP;
foreach $z (keys %branchhashP) {if ($branchhashP{$z}==1) {delete $branchhashP{$z}}}
$it=0;
while ($branchhashP>0)
{
$it++;
foreach $br (keys %branchhashP)
    {
    if (exists $branchhashPARS{$br}) 
	{if (exists $branchhash{$branchhashPARS{$br}}) 
	    {
	    print "$br -> $branchhashP{$br}"; 
	    $branchhashPARS{$branchhashP{$br}}=$branchhash{$branchhashPARS{$br}}; 
	    print "<<$branchhashPARS{$br}>> $branchhashP{$br} => $branchhash{$branchhashPARS{$br}}\n";
	    delete $branchhashP{$br}
	    }}
    elsif (exists $branchhash{$br} and exists $terminals{$br}) 
	{
	print "$br -> $branchhashP{$br}"; 
	$branchhashPARS{$branchhashP{$br}}=$branchhash{$br}; 
	print "<$branchhash{$br}> $branchhashP{$br} => $branchhash{$br}\n"; 
	delete $branchhashP{$br}
	}
    }
$branchhashP=keys %branchhashP;
print "\n\n\n$it $branchhashP ";
foreach $z (keys %branchhashP) {print "$z<$branchhashP{$z}> "}
print "\n\n";
}

%ancseqs=();
$ln=0;
open (PR, "outfile");
while (<PR>)
{
$ln++;
chomp;
if ($ln>6 and $_=~/\w/)
    {
    $_=~s/\s+$//;
    $_=~s/\?|\.|\*/X/g;
#    print STDERR "$_\n";
    while ($_=~/(\S+) (\S+)/) {$A=$1; $B=$2; $_=~s/$A $B/$A$B/;}
    $str=" $_";
    @str=split(/\s+/, $str);
#    print STDERR "$str[$#str]\n";
    if ($#str==4)
	{
	unless (exists $terminals{$str[2]})
	    {
	    if (exists $branchhashPARS{$str[2]})
		{
		$ancseqs{$branchhashPARS{$str[2]}}.=$str[$#str];
		}
	    }
	}
    else 
	{
	$ancseqs{ROOT}.=$str[$#str];
	}
    }
}
close PR;

open (O, ">$ARGV[2]");
foreach $anc (sort {$a cmp $b} keys %ancseqs)
{
print O "$anc $ancseqs{$anc}\n";
}
close O;
