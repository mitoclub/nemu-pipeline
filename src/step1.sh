#!/bin/bash

SPNAME=$1
SEQUENCE=$2
GENE=$3

DB=~/mtDB/mtDB
TBLASTN=~/bin/tblastn 

QUERY=query.fa
OUT=query.out-mitocode

printf ">$SPNAME\n$SEQUENCE\n" 1>$QUERY
if [ -e $QUERY ]
then
echo "$QUERY generated"
$TBLASTN -db $DB -db_gencode 2 -num_alignments 500 -query $QUERY -out $OUT 1>/dev/null 2>/dev/null
echo "tblastn run completed. $OUT generated"
fi
