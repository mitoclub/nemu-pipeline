#!/bin/bash

# inp
query_out_mitocode_sel_nuc=query.out-mitocode_sel.nuc

# out
query_out_mitocode_sel_unique_fa=query.out-mitocode_sel_unique.fa

script5=/opt/scripts/codon_alig_unique.pl

$script5 $query_out_mitocode_sel_nuc 1>$query_out_mitocode_sel_unique_fa 2>/dev/null
