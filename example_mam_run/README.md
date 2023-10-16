# Run pipeline

## NeMu pipeline inputs:

1. sequences (`mam_cytb_aln_sample.fasta`)
2. config file, that contains all input parameters (`comp-sp.config`)


## Execute pipeline:

```bash
nextflow -c comp-sp.config run nemu.nf --outdir outdir
```
