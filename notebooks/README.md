# Notebooks

1. [expectations_analysis.ipynb](./expectations_analysis.ipynb) - Expected mutations analysis. Contains (1) estimation of synonymous mutations counts in invariable positions; (2) calculation of expected mutations counts (used in calculation of mutation spectra) taking into account invariable positions
2. [mammals_spectra.ipynb](./mammals_spectra.ipynb) - Mammals spectra: compare different approaches of spectra reconstuction
3. [position_analysis.ipynb](./position_analysis.ipynb) - Selection exploration on Human and Mouse spectra (position-specific) on genes cytb, cox1, nd1 using novel method based on simulation of neutral evolution
4. [alisim_analysis_sp.ipynb](./alisim_analysis_sp.ipynb) - Analysis of simulated by `iqtree2 -alisim` random spectra on random trees (Figure 1 in the paper)
5. [alisim_analysis_mam.ipynb](./alisim_analysis_mam.ipynb) - Analysis of simulated by `iqtree2 -alisim` random spectra on big mammals tree (Figure 2 in the paper)
6. [prepare_beasts_datasets.ipynb](./prepare_beasts_datasets.ipynb) - Beast dataset preparation for NeMu pipeline execution. Extract mitochondrial genes from **MIDORI** database for different taxa and select outgroup sequence for each species-gene
7. [collect_spectra.ipynb](./collect_spectra.ipynb) - collect ann preanalyse spectra of mitochondrial genes of vertebrates
8. [prepare_nt_taxid.ipynb](./prepare_nt_taxid.ipynb) - analyse taxids of the nt database
9. [polymorphism_compare.ipynb](./polymorphism_compare.ipynb) - compare mitochondrial genes spectra of vertebrates, calculated by NeMu, with spectra, calculated by simple consensus-based approach
10. [prepare_taxa_table.ipynb](./prepare_taxa_table.ipynb) - - choose taxids that used as input for the pipeline on the webserver when nt database choosen
