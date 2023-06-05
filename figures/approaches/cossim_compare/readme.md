# Compare different spectra in *mammals*

GTR and RY10.12 are substitution models used in state reconstruction in [IQTREE2](http://www.iqtree.org/). 'Simple' suffix indicates spectra computed without states (and mutations) probabilities. [pastml](https://pastml.pasteur.fr/) label indicates spectra reconstderived with pastml tool (SYM model).


## Cutoffs search:

- [ratecat_cutoff_search_syn_pcutoff_0.1.pdf](ratecat_cutoff_search_syn_pcutoff_0.1.pdf) - positions of 3 and 4 rate-categories (medium substitution rate) too different from 5 and 6 rate-categories (high substitution rate)
- [proba_cutoff_search_syn_catcutoff_0.pdf](proba_cutoff_search_syn_catcutoff_0.pdf) - mutations with probaility < 0.3 are too different from other probability-categories. Used all sites.
- [proba_cutoff_search_syn_catcutoff_4.pdf](proba_cutoff_search_syn_catcutoff_4.pdf) - mutations with probaility < 0.3 are too different from other probability-categories. Used only 'fast' sites with high substitution rate.

## Compare different approaches of synonymous spectra adjusting (normalization) - COSMIC-like versus phylogenetically correct way:

- [pastml_compare_SYN_vs_SYN_C_pcutoff_0.3_catcutoff_4_mtypes16_mnum20_.pdf](pastml_compare_SYN_vs_SYN_C_pcutoff_0.3_catcutoff_4_mtypes16_mnum20_.pdf) - syn and syn_c are almost equal (99% of distribution have similarity > 0.8)

## Compare spectra derived using different methods of ancestral state recontruction with reference method **pastml**:

- [pastml_compare_pcutoff_0.0_catcutoff_0.pdf](pastml_compare_pcutoff_0.0_catcutoff_0.pdf) - Probability of mutations > 0 (P>0) and any sites 
- [pastml_compare_pcutoff_0.3_catcutoff_0.pdf](pastml_compare_pcutoff_0.3_catcutoff_0.pdf) - (P>0.3) and any sites
- [pastml_compare_pcutoff_0.3_catcutoff_4.pdf](pastml_compare_pcutoff_0.3_catcutoff_4.pdf) - (P>0.3) and only 'fast' sites (high substitution rate estimated in IQTREE2 using MLE)
- [pastml_compare_pcutoff_0.3_catcutoff_4_mtypes16_mnum20.pdf](pastml_compare_pcutoff_0.3_catcutoff_4_mtypes16_mnum20.pdf) - (P>0.3) and only 'fast' sites and filtration of tree edges with low number of mutations (<20) and low number of different mutation types (<16 out of 192)
