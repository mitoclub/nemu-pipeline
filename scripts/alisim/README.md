## Workflow

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
parallel --jobs 168 python ../../scripts/extract_inv_sites.py {} {.}.inv ::: generations_sp/*.fa

# extract real mutations
parallel --jobs 168 collect_mutations.py --tree external/mam_with_outgrp.nwk --states {} --states-fmt fasta --outdir spectra_groundtruth_mam/{/.} --gencode 2 --syn --force -q ::: generations_mam/*.fa
parallel --jobs 168 collect_mutations.py --tree {} --states {.}.fa --states-fmt fasta --outdir spectra_groundtruth_sp/{/.} --rates {.}.inv --gencode 2 --syn --force -q ::: generations_sp/*.nwk  # ???
# parallel --jobs 168 collect_mutations.py --tree {} --states {.}.fa --states-fmt fasta --outdir spectra_groundtruth_sp/{/.} --rates spectra_reconstructed_sp/{/.}/IQTREE/anc.rate --gencode 2 --syn --force -q ::: generations_sp/*.nwk

# run pipeline 
bash ./run_sp.sh
bash ./run_mam.sh

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
