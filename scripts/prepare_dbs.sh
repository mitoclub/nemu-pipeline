
cd ~/nemu-pipeline/data/MIDORI2/

# cd NUC
# parallel unzip {} ::: MIDORI2_UNIQ_NUC_*
# rm *.zip
# cd -

ls NUC/
# MIDORI2_UNIQ_NUC_GB259_A6_BLAST.fasta
# MIDORI2_UNIQ_NUC_GB259_A8_BLAST.fasta
# MIDORI2_UNIQ_NUC_GB259_CO1_BLAST.fasta
# MIDORI2_UNIQ_NUC_GB259_CO2_BLAST.fasta
# MIDORI2_UNIQ_NUC_GB259_CO3_BLAST.fasta
# MIDORI2_UNIQ_NUC_GB259_Cytb_BLAST.fasta
# MIDORI2_UNIQ_NUC_GB259_ND1_BLAST.fasta
# MIDORI2_UNIQ_NUC_GB259_ND2_BLAST.fasta
# MIDORI2_UNIQ_NUC_GB259_ND3_BLAST.fasta
# MIDORI2_UNIQ_NUC_GB259_ND4_BLAST.fasta
# MIDORI2_UNIQ_NUC_GB259_ND4L_BLAST.fasta
# MIDORI2_UNIQ_NUC_GB259_ND5_BLAST.fasta
# MIDORI2_UNIQ_NUC_GB259_ND6_BLAST.fasta

cd custom_db/
for file in ../NUC/MIDORI2_UNIQ_NUC_GB259_*; do
    interfasta=`basename $file`
    db_name=`basename $interfasta .fasta`
    echo $interfasta $db_name
    python3 ~/nemu-pipeline/scripts/process_midori_headers.py -i $file -o interim/$interfasta
    makeblastdb -in $interfasta -out $db_name -dbtype nucl -title "$db_name" -parse_seqids
done
cd -

# mv MIDORI2_UNIQ_NUC_* /scratch/blast_db/

