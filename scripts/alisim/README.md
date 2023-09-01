## Workflow

Here presented meta-pipeline of ALISIM analysis.

Stages:
1. Simulate trees and sequences including internal tree nodes genomes
2. Derive ground truth mutations and spectra from simulated data
3. Use terminal genomes from tree to reconstruct states and derive reconstructed mutations and spectra
4. Compare spectra from 2 and 3 

Bash code:

```bash
cd data/alisim

# add internal node names
scripts/assign_internal_nodes_labels.py data/external/timetrees/01_step2/4705sp_mean.nwk data/external/mam.nwk

# manually add OUTGRP node and Node4705 to the tree
# external/mam_with_outgrp.nwk

# generate alignments with internal genomes
scripts/ali_sim_run.py

# extract leaves from generated alignments
parallel --jobs 168 select_records.py -p Node -p ROOT -v {} {.}.fa.leaves ::: generations_mam/*.fa
parallel --jobs 168 select_records.py -p T -p OUTGRP  {} {.}.fa.leaves ::: generations_sp/*.fa

# extract inv sites from simulated sequences
parallel --jobs 168 python ../../scripts/alisim/extract_inv_sites.py {} {.}.inv ::: generations_sp/*.fa

# extract ground truth mutations
parallel --jobs 200 collect_mutations.py --tree external/mam_with_outgrp.nwk --states {} --states-fmt fasta --outdir spectra_groundtruth_mam/{/.} --gencode 2 --syn --force -q --save-exp-muts ::: generations_mam/*.fa
parallel --jobs 168 collect_mutations.py --tree {} --states {.}.fa --states-fmt fasta --outdir spectra_groundtruth_sp/{/.} --rates {.}.inv --gencode 2 --syn --force -q ::: generations_sp/*.nwk  # ???
# parallel --jobs 168 collect_mutations.py --tree {} --states {.}.fa --states-fmt fasta --outdir spectra_groundtruth_sp/{/.} --rates spectra_reconstructed_sp/{/.}/IQTREE/anc.rate --gencode 2 --syn --force -q ::: generations_sp/*.nwk

# run pipeline: reconstruct states and then spectra
bash ./run_sp.sh
bash ./run_mam.sh

# run spectra derivation after pipeline execution (need to compare correctly, so redo...)
cd spectra_reconstructed_mam
parallel --jobs 100 collect_mutations.py --tree {}/IQTREE/iqtree_anc_tree.nwk --states {}/IQTREE/iqtree_anc.state --states {}/tmp/leaves_states.state --gencode 2 --syn --proba --outdir {}/spectra_v3 --pcutoff 0.05 -q --save-exp-muts ::: *
# --rates {}/IQTREE/anc.rate --cat-cutoff 4
cd -

# delete tmp work files
for d in *; do rm -rf $d/work; done

```


## Run test:

```bash
cd test_sp_rates
nextflow -c  ../nemu_sp.config  run /home/kpotoh/nemu-pipeline/pipeline/Nemu-core/main.nf -with-trace --sequence ../generations_sp/gtr_100_rnd_replica_10.fa.leaves --outdir . --treefile ../generations_sp/gtr_100_rnd_replica_10.nwk
cd -
cd test_sp_ex_rates
nextflow -c  ../nemu_sp.config  run /home/kpotoh/nemu-pipeline/pipeline/Nemu-core/main.nf -with-trace --sequence ../generations_sp/gtr_100_rnd_replica_10.fa.leaves --outdir . --treefile ../generations_sp/gtr_100_rnd_replica_10.nwk --siterates.exclude_cons_sites false
```
