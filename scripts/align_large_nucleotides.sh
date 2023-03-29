#!/bin/bash

GENCODE=2
SUFFIX=".fna"
THREADS=12

FASTA=$1
PREFIX=$(basename $FASTA $SUFFIX)

# NT2AA
java -jar /home/kpotoh/tools/macse_v2.06.jar -prog translateNT2AA -seq $FASTA -gc_def $GENCODE -out_AA $PREFIX.faa
#ALN AA
mafft --thread $THREADS $PREFIX.faa > ${PREFIX}_aln.faa
#AA_ALN --> NT_ALN
java -jar /home/kpotoh/tools/macse_v2.06.jar -prog reportGapsAA2NT -align_AA ${PREFIX}_aln.faa -seq $FASTA -out_NT ${PREFIX}_aln.fna

rm $PREFIX.faa