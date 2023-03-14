## create nucleotide blast-db
makeblastdb -in humans_mt.fasta -dbtype nucl -title "Human mtDNA"
sudo singularity exec --bind /mnt/data/export:/export src/image_pipeline-2.5.sif makeblastdb -in /export/data/humans_mt.fasta -dbtype nucl -title "Human mtDNA" -parse_seqids


## Download refseq from 
https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/
