#!/bin/bash

# inp
SPNAME=$1
query_out_mitocode_sel_hash=query.out-mitocode_sel.hash

# out
query_out_mitocode_sel_nuc=query.out-mitocode_sel.nuc

script4=/opt/scripts/nuc_coding_mod.pl

printf "fasta protein alignment with $SPNAME and outgroup was generated\nand names legenda was printed (*.namehash)\n"
$script4 $query_out_mitocode_sel_hash 1>$query_out_mitocode_sel_nuc 2>/dev/null
