## Scripts used in NeMu-pipeline

- [codon_alig_unique.pl](./codon_alig_unique.pl) - drop duplicated sequences
- [header_sel_mod3.pl](./header_sel_mod3.pl) - selection of species and outgroup sequences for analysis
- [macse2.pl](./macse2.pl) - QC after macse alignment
- [mut_calculation_step.sh](./mut_calculation_step.sh) - legacy??? part of mut extraction step 
- [mutnumbers.pl](./mutnumbers.pl) - estimating maximum theoretical number of variable positions and mutations
- [nuc_coding_mod.pl](./nuc_coding_mod.pl) - nucleotide sequences extraction from BLAST database
- [rename_tree_nodes.sh](./rename_tree_nodes.sh) - part of simulation (pyvolve) process 
- [resci.py](./resci.py) - reformat scientific notation to float 
- [blasting.sh](./blasting.sh) - tblastn for homologous sequences search 

## Scripts-formatters

- [aln2pastml.py](./aln2pastml.py) - prepare dataset with several genes for pastml run
- [aln2pastml_single.py](./aln2pastml_single.py) - prepare dataset with single gene for pastml run
- [fasta2fasta_2line.py](./fasta2fasta_2line.py) - reformat fasta multiline to 2line style
- [gbk2fasta.py](./gbk2fasta.py) - full sequence from genbank to fasta
- [pastml2custom_format.py](./pastml2custom_format.py) - pastml output to iqtree states format with Part column for several genes
- [pastml2custom_format_single.py](./pastml2custom_format_single.py) - pastml output to iqtree states format with Part column for single gene

## Processing

- [align_large_nucleotides.sh](./align_large_nucleotides.sh) - 
- [exposure_completing.sh](./exposure_completing.sh) - selection extraction using evolution simulation
- [prepare_ND1_rate.py](./prepare_ND1_rate.py) - 
- [process_midori_headers.py](./process_midori_headers.py) - extract species names from MIDORI records and place them after record index (>ID NAME DESCR)
- [run_pyvolve.sh](./run_pyvolve.sh) - legacy??? TODO 
- [split_genomes_to_genes.py](./split_genomes_to_genes.py) - split GAGP dataset to genes
- [ali_sim_run.py](./ali_sim_run.py) - Run simulation to get test data for pipeline verification 


## External scripts from [PyMutSpec](https://github.com/mitoclub/PyMutSpec)

Many [scripts](https://github.com/mitoclub/PyMutSpec/tree/master/scripts) and [functionalities](https://github.com/mitoclub/PyMutSpec/tree/master/pymutspec) in NeMu-pipeline and data analysis stored in separated Python library - PyMutSpec
