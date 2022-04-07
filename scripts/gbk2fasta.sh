#!/bin/bash

for file in $@ ; do
    echo $file
    python3.8 -c "from Bio import SeqIO; SeqIO.convert(\"$file\", genbank, ${file%.*}.fasta, fasta);"
done
