#!/bin/bash

# inp
query_out_fasta=query.out-mitocode.fa

# out
query_out_mitocode_sel_fa=query.out-mitocode_sel.fa
query_out_mitocode_sel_hash=query.out-mitocode_sel.hash

script3=/opt/scripts/header_sel_mod3.pl

$script3 $query_out_fasta $SPNAME 1>$query_out_mitocode_sel_fa 2>$query_out_mitocode_sel_hash

