#!/bin/bash

THREADS=16
gencode=2
max_target_seqs=1500
required_nseqs=3
outfmt="6 saccver pident length qlen gapopen sstart send evalue bitscore sframe"
DB=/home/dolphin/db/nt

species_name='Mus musculus'
genus_taxid=10088
query=mouse_cytb.faa


echo "INFO: Collecting genus taxids" >&2
if [[ "$species_name" == Homo* ]]; then
    get_species_taxids.sh -t 9604 > genus.taxids
else
    get_species_taxids.sh -t $genus_taxid > genus.taxids
fi

sleep 2
echo "INFO: Collecting species taxid information" >&2
get_species_taxids.sh -n "$species_name" | head -n 5 > sp_tax.info
if [[ `grep "rank : species" sp_tax.info` ]]; then 
    raw_sp_taxid=`grep Taxid sp_tax.info`
    species_taxid="${raw_sp_taxid#*Taxid : }"
fi
sleep 1
echo "INFO: Collecting under-species taxids" >&2
get_species_taxids.sh -t $species_taxid > species.taxids
# echo $species_taxid > species.taxids

grep -v -f species.taxids genus.taxids > relatives.taxids

echo "INFO: Checking number of taxids for outgroup" >&2
if [ `wc -l relatives.taxids | cut -f 1 -d ' '` -eq 0 ]; then
    echo "ERROR: there are no taxids that can be used as outgroup." >&2
    echo "Maybe this species is single in the genus, so pipeline can build incorrect phylogeneti tree" >&2
    echo "due to potential incorrect tree rooting. You can select sequences and outgroup manually and" >&2
    echo "run pipeline on your nucleotide sequences" >&2
    exit 1
fi

echo "INFO: Blasting species sequences in the nt" >&2
tblastn -db $DB -db_gencode $gencode -max_target_seqs $max_target_seqs \
        -query $query -out blast_output_species.tsv -evalue 0.00001 \
        -num_threads $THREADS -taxidlist species.taxids \
        -outfmt "$outfmt"

echo "INFO: Filtering out bad hits: ident <= 70, query coverage <= 0.6" >&2
awk '$2 > 70 && $3 / $4 > 0.6' blast_output_species.tsv > blast_output_species_filtered.tsv

echo "INFO: Checking required number of hits" >&2
nhits=`wc -l blast_output_species_filtered.tsv | cut -f 1 -d ' '`
if [ $nhits -lt $required_nseqs ]; then
    echo "ERROR: there are only $nhits valuable hits in the database for given query," >&2
    echo "but needed at least $required_nseqs" >&2
    exit 1
fi

echo "INFO: Preparing coords for nucleotide sequences extraction" >&2
# entry|range|strand, e.g. ID09 x-y minus
awk '{print $1, $6, $7, ($NF ~ /^-/) ? "minus" : "plus"}' blast_output_species_filtered.tsv > raw_coords.txt
awk '$2 > $3 {print $1, $3 "-" $2, $4}' raw_coords.txt > coords.txt
awk '$3 > $2 {print $1, $2 "-" $3, $4}' raw_coords.txt >> coords.txt

echo "INFO: Getting nucleotide sequences" >&2
blastdbcmd -db $DB -entry_batch coords.txt -outfmt %f -out nucleotide_sequences.fasta

echo "INFO: Checking required number of extracted seqs" >&2
nseqs=`grep -c '>' nucleotide_sequences.fasta`
if [ $nseqs -lt $required_nseqs ]; then
    echo "ERROR: cannot extract more than $nseqs seqs from the database for given query, but needed at least $required_nseqs" >&2
    exit 1
fi

echo -e "INFO: Blasting for outgroup search" >&2
tblastn -db $DB -db_gencode $gencode -max_target_seqs 10 \
        -query $query -out blast_output_genus.tsv -evalue 0.00001 \
        -num_threads $THREADS -taxidlist relatives.taxids \
        -outfmt "$outfmt"

echo "INFO: Filtering out bad hits: ident <= 70, query coverage <= 0.6" >&2
awk '$2 > 70 && $3 / $4 > 0.6' blast_output_genus.tsv | sort -rk 9 > blast_output_genus_filtered.tsv

echo "INFO: Checking required number of hits for outgroup" >&2
if [ `wc -l blast_output_genus_filtered.tsv | cut -f 1 -d ' '` -eq 0 ]; then
    echo "ERROR: there are no hits in the database that could be used as outgroup" >&2
    exit 1
fi

echo "INFO: Preparing genus coords for nucleotide sequences extraction" >&2
# entry|range|strand, e.g. ID09 x-y minus
head -n 1 blast_output_genus_filtered.tsv | awk '{print $1, $6, $7, ($NF ~ /^-/) ? "minus" : "plus"}' > raw_coords_genus.txt
awk '$2 > $3 {print $1, $3 "-" $2, $4}' raw_coords_genus.txt >  coords_genus.txt
awk '$3 > $2 {print $1, $2 "-" $3, $4}' raw_coords_genus.txt >> coords_genus.txt

echo "INFO: Getting nucleotide sequences of potential outgroups" >&2
blastdbcmd -db $DB -entry_batch coords_genus.txt -outfmt %f -out outgroup_sequence.fasta

echo "INFO: Sequences headers encoding" >&2
ohead=`head -n 1 outgroup_sequence.fasta`
cat outgroup_sequence.fasta nucleotide_sequences.fasta > sample.fasta
