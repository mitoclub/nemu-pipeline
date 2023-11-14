
THREADS=4
gencode=2

# SPECIES='Stenella longirostris'
SPECIES='Homo sapiens'

gene=Cytb
DB=/scratch/blast_db/MIDORI2_UNIQ_NUC_GB253_${gene}_BLAST
query=$1

report=report.blast
nseqs=500
LIMIT=33000

# OUTGROUP QC
while true
do   
    echo "INFO: run blasting; nseqs=${nseqs}"
    tblastn -db $DB -db_gencode $gencode -num_descriptions $nseqs -num_alignments $nseqs -query $query -out $report -num_threads $THREADS
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
            echo "UNEXPECTED ERROR: database cannot contain more than $LIMIT sequences of one species"
            exit 1
        fi
        echo "INFO: run blasting again due to abcence of different species; nseqs=${nseqs}"
    else
        echo "SUCCESS: other species for outgroup selection is in hits"
        break
    fi
done
