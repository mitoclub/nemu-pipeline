
THREADS=1
gencode=2

# SPECIES='Homo sapiens'
SPECIES='Stenella longirostris'
GENUS_TAXID=9734

DB=/home/dolphin/db/nt
query=/home/kpotoh/nemu-pipeline/example_data/CYTB__Stenella_longirostris.faa

report=report.blast
nseqs=500
LIMIT=33000

get_species_taxids.sh -t $GENUS_TAXID > genus_species_taxids.txt

TAXIDS=`cat genus_species_taxids.txt | paste -sd,`
echo -e "Search in these taxids:\n$TAXIDS\n"

while true
do   
    echo "INFO: run blasting; nseqs=${nseqs}"
        tblastn -db $DB -db_gencode $gencode -num_descriptions $nseqs -num_alignments $nseqs \
                -query $query -out $report -evalue 0.001 -num_threads $THREADS -taxids $TAXIDS

    if [ `grep -c "No hits found" $report` -eq 0 ]; then 
        echo "SUCCESS: hits found in the database for given query"
    else
        echo "ERROR: there are no hits in the database for given query" > no_hits.log
        cat no_hits.log
        exit 1
    fi

    if [ `grep -c "$SPECIES" $report` -ge $((nseqs * 2 - 10)) ]; then
        nseqs=$((nseqs * 4))
        if [ $nseqs -gt $LIMIT ]; then
            echo "UNEXPECTED ERROR: database cannot contain more than $LIMIT sequences of one gene of species"
            exit 1
        fi
        echo "INFO: run blasting again due to abcence of different species; nseqs=${nseqs}"
    else
        echo "SUCCESS: other species for outgroup selection is in hits"
        break
    fi
done
