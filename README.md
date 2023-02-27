# Nemu-pipeline

The pipeline for neutral mutation spectra evaluation based on evolutionary data.

We use [Dolphinnext]() system for interctive user-friendly execution of the pipeline.

https://biopipelines.kantiana.ru/dolphinnext


## create nucleotide blast-db
makeblastdb -in humans_mt.fasta -dbtype nucl -title "Human mtDNA"
sudo singularity exec --bind /mnt/data/export:/export src/image_pipeline-2.5.sif makeblastdb -in /export/data/humans_mt.fasta -dbtype nucl -title "Human mtDNA" -parse_seqids


## containers
1. p.def - only python and its programms
2. pipeline-2.5.def - working through p.def result container
3. pipeline.def - independent and working


## Download refseq from 
https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/
