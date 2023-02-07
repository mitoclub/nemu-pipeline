#!/usr/bin/perl -w

%seqs=();
@seqs=();
open (I, "$ARGV[0]");
while (<I>)
{
chomp;
#fasta header
if ($_=~/^>(.*?)$/)
    {
    $name=$1;
    push @seqs, $name;
    }
#fasta sequence
else
    {
    @I=split(//, $_);
    $codonPos=0;
    $codon="";
    @codons=();
    for ($j=0;$j<$#I+1;$j++)
	{
	$codonPos++;
	$codon.=$I[$j];
	if ($codonPos==3)
	    {
	    push @codons, $codon;
	    $codonPos=0;
	    $codon="";
	    }
	}
    @{$seqs{$name}}=@codons;
    $max=$#codons+1;
    }
}
close I;

%nonnuccodons=();
for ($i=0;$i<$max;$i++)
{
$nonnuccodons{$i}=0;
foreach $nm (keys %seqs){
    if ($seqs{$nm}[$i]=~/[^ACGTacgt]/) {
		print $seqs{$nm}[$i];
		$nonnuccodons{$i}+=1;
		if ($seqs{$nm}[$i]=~/[!\?\.,]/) {
			$seqs{$nm}[$i]="NNN"
		}
	}
}
}

#if 30% of codons in alignment column than remove this column
$nseqs = @seqs;
$thresh=int(0.3 * $nseqs);
print "nseqs = $nseqs, threshold = $thresh\n";

open (O, ">$ARGV[1]");
foreach $seq (@seqs)
{
print O ">$seq\n";
for ($i=0;$i<$max;$i++) {
	if ($nonnuccodons{$i} < $thresh) {print O "$seqs{$seq}[$i]"}
}
print O "\n";
}
close O;
