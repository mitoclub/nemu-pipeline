# *Mammalia* workflow

1. Build tree

...

2. Reconstruct ancestor sequences

pastml workflow [here](./pastml_readme.md)

```bash
# TODO add pastml commands

nohup iqtree2 -te iqtree_anc_tree.nwk -s alignment_checked.fasta -m GTR+FO+G6+I -asr -nt 64 --prefix states/gtr --rate &
nohup iqtree2 -te iqtree_anc_tree.nwk -s alignment_checked.fasta -m STRSYM+FO+G6+I -asr -nt 64 --prefix states/strsym --rate &
nohup iqtree2 -te iqtree_anc_tree.nwk -s alignment_checked.fasta -m UNREST+FO+G6+I -asr -nt 64 --prefix states/unrest --rate &

# and slightly reformat
iqtree_states_add_part.py gtr.state gtr.state.custom
iqtree_states_add_part.py unrest.state unrest.state.custom
iqtree_states_add_part.py strsym.state strsym.state.custom
```

3. Extract mutations from trees (TODO many descriptrion... PyMutSpec + simple + proba...)

```bash
#pastml states
nohup collect_mutations.py --tree iqtree_anc_tree.nwk --states ../pastml_states/mammals_cytb_pastml_states.tsv --rates states/unrest.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --proba --no-phylocoef --outdir pastml &
nohup collect_mutations.py --tree iqtree_anc_tree.nwk --states ../pastml_states/mammals_nd1_pastml_states.tsv --rates states/unrest.rate  --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --proba --no-phylocoef --outdir pastml &

#iqtree states
collect_mutations.py --tree iqtree_anc_tree.nwk --states states/gtr.state.custom --states leaves_states.state --rates states/gtr.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --outdir gtr_simple
collect_mutations.py --tree iqtree_anc_tree.nwk --states states/gtr.state.custom --states leaves_states.state --rates states/gtr.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --proba --outdir gtr_proba
collect_mutations.py --tree iqtree_anc_tree.nwk --states states/strsym.state.custom --states leaves_states.state --rates states/strsym.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --outdir strsym_simple
collect_mutations.py --tree iqtree_anc_tree.nwk --states states/strsym.state.custom --states leaves_states.state --rates states/strsym.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --proba --outdir strsym_proba
collect_mutations.py --tree iqtree_anc_tree.nwk --states states/unrest.state.custom --states leaves_states.state --rates states/unrest.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --outdir unrest_simple
collect_mutations.py --tree iqtree_anc_tree.nwk --states states/unrest.state.custom --states leaves_states.state --rates states/unrest.rate --gencode 2 --no-mutspec --save-exp-muts --syn --syn_c --syn4f --proba --outdir unrest_proba
```
