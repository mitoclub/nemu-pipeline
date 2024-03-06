# Nemu-pipeline

The pipeline for neutral mutation spectra evaluation based on evolutionary data.

Check [Wiki](https://github.com/mitoclub/nemu-pipeline/wiki) to read how to use the pipeline

## Webserver

Enjoy the pipeline using our user-friendly [webserver](https://nemu-pipeline.streamlit.app/)!

## Nextflow script

NeMu pipeline developed using Nextflow. Available [here](./pipeline/).

## Dependencies

## Singularity container

The pipeline relies on a diverse set of dependencies. To simplify execution, ensure compatibility across different operating systems, all these dependencies have been packaged within a [Singularity container](./singularity/).

### Used software

Here listed all languages, programs and packages used in the NeMu pipeline during mutatoinal spectra calculations:

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

## Article

Read our draft article for more details!

- NeMu: A Comprehensive Pipeline for Accurate Reconstruction of Neutral Mutation Spectra 
from Evolutionary Data Bogdan Efimenko, Konstantin Popadin, Konstantin Gunbin 
**bioRxiv** 2023.12.13.571433; doi: https://doi.org/10.1101/2023.12.13.571433
