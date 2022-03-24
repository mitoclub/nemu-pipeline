#!/bin/bash

# inp
SPNAME=$1
SEQUENCE=$2
GENE=$3

DB=~/mtDB/mtDB
TBLASTN=~/bin/tblastn 

query=query.fa

# out
query_out=query.out-mitocode

printf ">$SPNAME\n$SEQUENCE\n" 1>$query

echo "$query generated"
$TBLASTN -db $DB -db_gencode 2 -num_alignments 500 -query $query -out $query_out 1>/dev/null 2>/dev/null
# echo "tblastn run completed. $query_out generated"
# rm $query
