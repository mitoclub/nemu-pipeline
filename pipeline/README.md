# The pipeline for neutral mutation spectra evaluation based on evolutionary data

## 2 pipeline versions

1. [NeMu-pipeline including tblastn-head](./nemu.nf) - input is single protein sequence that will be used by tblastn to search homologous nucleotide sequences in selected database. After this step NeMu-core executes with phylogeny and spectra inference
2. [NeMu-core pipeline](./nemu-core.nf) - input is multifasta of nucleotide sequences, that used for phylogeny and spectra inference

## Config examples

1. [Config for comparative species analysis](./comp_sp.config) - on many species
2. [Config for intraspecies analysis](./single_sp.config) - on single species

    - Don't forget to change process.container and singularity.runOptions parameters according to execution environment
