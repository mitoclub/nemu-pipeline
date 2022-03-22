#!/bin/bash

# inp
query_out=query.out-mitocode

# out
query_out_fasta=query.out-mitocode.fa

MVIEW=mview

echo "protein homologs were found"
$MVIEW -in blast -out fasta $query_out 1>$query_out_fasta 2>/dev/null
