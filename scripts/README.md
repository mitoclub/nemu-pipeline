## Scripts used in NeMu pipeline

- [blasting.sh](./blasting.sh) - tblastn for homologous sequences search in midori databases (integrated directly to the pipeline script)
- [blasting_nt.sh](./blasting_nt.sh) - tblastn for homologous sequences search in nt database (integrated directly to the pipeline script)
- [perl/nuc_coding_mod.pl](./perl/nuc_coding_mod.pl) - nucleotide sequences extraction from BLAST database
- [perl/macse2.pl](./perl/macse2.pl) - QC after macse alignment
- [perl/codon_alig_unique.pl](./perl/codon_alig_unique.pl) - drop duplicated sequences
- [perl/header_sel_mod3.pl](./perl/header_sel_mod3.pl) - selection of species and outgroup sequences for analysis
- [legacy/resci.py](./legacy/resci.py) - reformat scientific notation to float (legacy)
- [mut_calculation_step.sh](./mut_calculation_step.sh) - legacy, part of mut extraction step
- [rename_tree_nodes.sh](./rename_tree_nodes.sh) - part of simulation (pyvolve) process

## Scripts-formatters

- [aln2pastml.py](./aln2pastml.py) - prepare dataset with several genes for pastml run
- [aln2pastml_single.py](./aln2pastml_single.py) - prepare dataset with single gene for pastml run
- [fasta2fasta_2line.py](./fasta2fasta_2line.py) - reformat fasta multiline to 2line style
- [gbk2fasta.py](./gbk2fasta.py) - full sequence from genbank to fasta
- [pastml2custom_format.py](./pastml2custom_format.py) - format pastml output to iqtree states format with Part column for several genes
- [pastml2custom_format_single.py](./pastml2custom_format_single.py) - pastml output to iqtree states format with Part column for single gene

## Processing

- [ali_sim_run.py](./ali_sim_run.py) - Run simulation to get test data for pipeline verification 
- [run_simmulation.sh](./run_simmulation.sh) - simulate neutral evolution on human and murine COX1, Cytb and ND1 genes  
- [align_large_nucleotides.sh](./align_large_nucleotides.sh) - 
- [prepare_ND1_rate.py](./prepare_ND1_rate.py) - 
- [process_midori_headers.py](./process_midori_headers.py) - extract species names from MIDORI records and place them after record index (>ID NAME DESCR)
- [split_genomes_to_genes.py](./split_genomes_to_genes.py) - split GAGP dataset to genes


## External scripts from [PyMutSpec](https://github.com/mitoclub/PyMutSpec)

Many [scripts](https://github.com/mitoclub/PyMutSpec/tree/master/scripts) and [functionalities](https://github.com/mitoclub/PyMutSpec/tree/master/pymutspec) in NeMu pipeline and data analysis stored in separated Python library - PyMutSpec
