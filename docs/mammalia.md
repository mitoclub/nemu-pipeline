# *Mammalia* workflow

1. Build tree

...

2. Reconstruct ancestor sequences

pastml workflow [here](./pastml_readme.md)

```bash
# TODO add pastml commands

iqtree2 -te iqtree_anc_tree.nwk -s alignment_checked.fasta -m GTR+FO+G6+I -asr -nt 64 --prefix states/gtr --rate
iqtree2 -te iqtree_anc_tree.nwk -s alignment_checked.fasta -m STRSYM+FO+G6+I -asr -nt 64 --prefix states/strsym --rate
# iqtree2 -te iqtree_anc_tree.nwk -s alignment_checked.fasta -m UNREST+FO+G6+I -asr -nt 64 --prefix states/unrest --rate # ERROR SPECTRA
iqtree2 -te iqtree_anc_tree.nwk -s alignment_checked.fasta -m 6.8a+FO+G6+I -asr -nt 64 --prefix states/RY6.8a --rate
iqtree2 -te iqtree_anc_tree.nwk -s alignment_checked.fasta -m 8.8+FO+G6+I -asr -nt 64 --prefix states/RY8.8 --rate
iqtree2 -te iqtree_anc_tree.nwk -s alignment_checked.fasta -m 10.12+FO+G6+I -asr -nt 64 --prefix states/RY10.12 --rate

# and slightly reformat
iqtree_states_add_part.py gtr.state gtr.state.custom
iqtree_states_add_part.py unrest.state unrest.state.custom
iqtree_states_add_part.py strsym.state strsym.state.custom
iqtree_states_add_part.py RY6.8a.state RY6.8a.state.custom
iqtree_states_add_part.py RY8.8.state RY8.8.state.custom
iqtree_states_add_part.py RY10.12.state RY10.12.state.custom
```

3. Extract mutations from trees (TODO many descriptrion... PyMutSpec + simple + proba...)

```bash
#pastml states
collect_mutations.py --tree iqtree_anc_tree.nwk --states ../pastml_states/mammals_cytb_pastml_states.tsv --rates states/unrest.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --proba --no-phylocoef --outdir pastml
collect_mutations.py --tree iqtree_anc_tree.nwk --states ../pastml_states/mammals_nd1_pastml_states.tsv --rates states/unrest.rate  --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --proba --no-phylocoef --outdir pastml

#iqtree states
for model in gtr strsym unrest RY6.8a RY8.8 ; do
    collect_mutations.py --tree iqtree_anc_tree.nwk --states states/${model}.state.custom --states leaves_states.state --rates states/${model}.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --outdir ${model}_simple
    collect_mutations.py --tree iqtree_anc_tree.nwk --states states/${model}.state.custom --states leaves_states.state --rates states/${model}.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --outdir ${model}_proba --proba
done



collect_mutations.py --tree iqtree_anc_tree.nwk --states states/gtr.state.custom --states leaves_states.state --rates states/gtr.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --outdir gtr_simple
collect_mutations.py --tree iqtree_anc_tree.nwk --states states/gtr.state.custom --states leaves_states.state --rates states/gtr.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --proba --outdir gtr_proba
collect_mutations.py --tree iqtree_anc_tree.nwk --states states/strsym.state.custom --states leaves_states.state --rates states/strsym.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --outdir strsym_simple
collect_mutations.py --tree iqtree_anc_tree.nwk --states states/strsym.state.custom --states leaves_states.state --rates states/strsym.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --proba --outdir strsym_proba
collect_mutations.py --tree iqtree_anc_tree.nwk --states states/unrest.state.custom --states leaves_states.state --rates states/unrest.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --outdir unrest_simple
collect_mutations.py --tree iqtree_anc_tree.nwk --states states/unrest.state.custom --states leaves_states.state --rates states/unrest.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --proba --outdir unrest_proba

collect_mutations.py --tree iqtree_anc_tree.nwk --states states/RY6.8a.state.custom --states leaves_states.state --rates states/RY6.8a.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --outdir RY6.8a_simple
collect_mutations.py --tree iqtree_anc_tree.nwk --states states/RY6.8a.state.custom --states leaves_states.state --rates states/RY6.8a.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --proba --outdir RY6.8a_proba

collect_mutations.py --tree iqtree_anc_tree.nwk --states states/RY8.8.state.custom --states leaves_states.state --rates states/RY8.8.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --outdir RY8.8_simple
collect_mutations.py --tree iqtree_anc_tree.nwk --states states/RY8.8.state.custom --states leaves_states.state --rates states/RY8.8.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --proba --outdir RY8.8_proba

nohup collect_mutations.py --tree iqtree_anc_tree.nwk --states states/RY10.12.state.custom --states leaves_states.state --rates states/RY10.12.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --outdir RY10.12_simple &
nohup collect_mutations.py --tree iqtree_anc_tree.nwk --states states/RY10.12.state.custom --states leaves_states.state --rates states/RY10.12.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --proba --outdir RY10.12_proba &

```
