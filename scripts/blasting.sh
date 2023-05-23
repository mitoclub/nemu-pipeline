
THREADS=16
gencode=2
nseqs=1000

gene=Cytb
SPECIES="Homo sapiens"

DB=/scratch/blast_db/MIDORI2_UNIQ_NUC_GB253_${gene}_BLAST
query=$1

report=report.blast

tblastn -db $DB -db_gencode $gencode -num_descriptions $nseqs -num_alignments $nseqs -query $query -out $report -num_threads $THREADS

if [ `grep -c "No hits found" report.blast` -eq 0 ]; then 
	echo "SUCCESS: hits found in the database for given query"
else
	echo "ERROR: there are no hits in the database for given query" > no_hits.log
	cat no_hits.log
	exit 1
fi

# OUTGROUP QC

if [ `grep -c $SPECIES $report` -ge (nseqs * 2)TODO ]; then
    ...
    nseqs=5000
    run blasting again
    write it like for loop with 'exit 0'

fi