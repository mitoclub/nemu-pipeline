# Nemu-pipeline

[![NeMu Paper](https://img.shields.io/badge/DOI-10.1093%2Fnar%2Fgkae438-blue)](https://doi.org/10.1093/nar/gkae438)
[![Static Badge](https://img.shields.io/badge/Nextflow-lightblue)](https://www.nextflow.io/)
[![Static Badge](https://img.shields.io/badge/Singularity-green)](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html)
![GitHub commit activity](https://img.shields.io/github/commit-activity/y/mitoclub/nemu-pipeline)
![GitHub last commit](https://img.shields.io/github/last-commit/mitoclub/nemu-pipeline)


This repository contains the pipeline for neutral mutational spectra evaluation based on evolutionary data and materials for the publication titled "NeMu: A Comprehensive Pipeline for Accurate Reconstruction of Neutral Mutation Spectra from Evolutionary Data" by Efimenko B. et al.

## Webserver

Enjoy the pipeline using our user-friendly [webserver](https://nemu-pipeline.com/)!

Check [Wiki](https://nemu-pipeline.com/Information) to read how to use the pipeline

## Repository details

- `./pipeline/` directory contains *.nf files (Nextflow) with 2 versions of the pipeline: for single-protein and multiple-nucleotide inputs
- `./singularity/` - definition files for the container used to store all dependencies
- `./data/` - several intermediate files from analyses, example inputs, full example run and supplementary data for the article
- `./figures/` - figures from all performed analyses and comparisons, pipeline schemes, etc.
- `./notebooks/` - all performed analyses on python
- `./scripts/` - all used scripts for specific data processing during analyses and scripts, used in the NeMu pipeline
- `./docs/` - notes for used tools

Each main directory contains readme file with detailed information.

## Dependencies

### Singularity container

The pipeline relies on a diverse set of dependencies. To simplify execution, ensure compatibility across different operating systems, all these dependencies have been packaged within a [Singularity container](./singularity/). We used Singularity 3.10.2.

### Used software

Here listed all dependencies used in the NeMu pipeline during mutatoinal spectra calculations:

- Linux (e.g. Ubuntu)
- Nextflow 22.10.7
- Python 3.8.12
- Perl v5.16.3
- R 3.6.0

1. **BLAST+** 2.13.0 for tblast search in nucleotide sequences databases;
2. **Mview** 1.67 for blast-report reformatting;
3. **Custom Perl scripts** for sequences extraction and filtration;
4. **MACSE** v2.06 and **MAFFT** 7.505 for nucleotide sequences alignment considering codon structure;
5. **IQ-TREE2** 2.2.0 for tree building, ancestral reconstruction and site rates estimation;
6. **TreeShrink** 1.3.9 for outlier branch filtration;
7. **Newick-tools** 1.6 for tree rerooting and other tree processing;
8. **Pymutspec** 0.0.8 Python package that we developed for mutation extraction and spectrum calculation.
9. **Pyvolve** 1.1.0 (modified) for neutral evolution simutation using MutSel models

## How to cite?

```text
Bogdan Efimenko, Konstantin Popadin, Konstantin Gunbin, NeMu: a comprehensive pipeline for accurate reconstruction of neutral mutation spectra from evolutionary data, Nucleic Acids Research, Volume 52, Issue W1, 5 July 2024, Pages W108–W115, https://doi.org/10.1093/nar/gkae438
```