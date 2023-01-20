
indir=...

mkdir $indir/ms && python3 ~/mutspec-utils/scripts/calculate_mutspec.py -b ${indir}/observed_mutations_iqtree.tsv -e ${indir}/exp_muts_invariant.tsv -o $indir/ms \
    --exclude OUTGRP,ROOT --proba  --syn --plot -x pdf


