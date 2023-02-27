
## Extration of mutations from tree

```bash
collect_mutations.py --tree iqtree_anc_tree.nwk --states iqtree_anc.state --states leaves_states.state --gencode 2 --syn --syn4f --outdir simple
collect_mutations.py --tree iqtree_anc_tree.nwk --states iqtree_anc.state --states leaves_states.state --gencode 2 --syn --syn4f --proba --outdir proba
collect_mutations.py --tree iqtree_anc_tree.nwk --states iqtree_anc.state --states leaves_states.state --gencode 2 --syn --syn4f --proba --no-phylocoef --outdir proba-no-phylocoef
# pastml mut extraction [here](../data/exposure/pastml_states/process_states.sh)
```