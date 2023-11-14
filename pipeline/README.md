# The pipeline for neutral mutation spectra evaluation based on evolutionary data

## How to run it on custom dataset

### Prepare dataset

1. Each gene of each species must be in separate fasta file
2. Each fasta file must contain one outgroup sequence (most closest gene sequence from another species) with header ">OUTGRP"
3. Fasta files must have standard extension ".fasta". ".fna", ".fa" and other forms are not appropriate
4. Files must be in one directory

### Run pipeline on prepared dataset

1. Copy [pipeline file](../Nemutail/main.nf) and [config file](../Nemutail/nextflow.config) to dataset directory
2. Modify config file according to your computer resources
3. Run `parallel nextflow -q -log {\.}.log run main.nf --sequence {} ::: *.fasta` inside dataset directory
4. Pipeline will create directories for each run

Dataset directory after execution will looks like:
```
$ cd dataset/
$ nextflow main.nf --sequence *.fasta
$ ls
A6__Dascyllus_aruanus
A6__Dascyllus_reticulatus
A6__Dawkinsia_denisonii
A6__Delphinapterus_leucas
A6__Dendrocoptes_medius
A6__Dascyllus_aruanus.fasta
A6__Dascyllus_reticulatus.fasta
A6__Dawkinsia_denisonii.fasta
A6__Delphinapterus_leucas.fasta
A6__Dendrocoptes_medius.fasta
main.nf
nextflow.config
work
```
