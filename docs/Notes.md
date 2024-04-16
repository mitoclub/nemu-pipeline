## create nucleotide blast-db
makeblastdb -in humans_mt.fasta -dbtype nucl -title "Human mtDNA"
sudo singularity exec --bind /mnt/data/export:/export src/image_pipeline-2.5.sif makeblastdb -in /export/data/humans_mt.fasta -dbtype nucl -title "Human mtDNA" -parse_seqids


## Download refseq from 
https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/


### Installing newick-utils-1.6

```bash
wget http://bioinfodbs.kantiana.ru/newick-utils-1.6.tar.gz
tar -xvzf newick-utils-1.6.tar.gz
cd newick-utils-1.6
./configure --prefix=/opt/newick-utils-1.6
make install
cd build/bin
ls
```