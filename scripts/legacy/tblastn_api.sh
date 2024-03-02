#!/bin/bash

# Check for proper number of command line args.
EXPECTED_ARGS=3
E_BADARGS=65

if [ $# -ne $EXPECTED_ARGS ]; then
    echo "Usage: $0 fasta_file genetic_code_number species_name"
    exit $E_BADARGS
fi

# Assign the arguments to variables
FASTA_FILE=$1
GENETIC_CODE=$2
SPECIES_NAME=$3

# Extract the protein sequence from the FASTA file
PROTEIN_SEQUENCE=$(awk '/^>/ {if (seqlen){print seqlen}; print;seqlen=0;next; } { printf "%s", $0;seqlen+=length($0);} END {print seqlen;}' $FASTA_FILE | tail -n +2)

# Send the query to the NCBI blast api
# Here we use tblastn, which translates a protein query sequence and compares it to a nucleotide sequence database dynamically translated in all reading frames.
RESPONSE=$(curl -s -H "Content-Type: application/x-www-form-urlencoded" -d "CMD=Put&PROGRAM=tblastn&DATABASE=nr&QUERY=${PROTEIN_SEQUENCE}&GENETIC_CODE=${GENETIC_CODE}&ORGANISM=${SPECIES_NAME}" https://blast.ncbi.nlm.nih.gov/Blast.cgi)

# Extract the Request ID (RID) and Estimated Time (RTOE) from the response
RID=$(echo $RESPONSE | grep -oP 'RID = \K\w+')
RTOE=$(echo $RESPONSE | grep -oP 'RTOE = \K\d+')

# Wait for search to complete before trying to get the results
# NCBI recommends waiting RTOE * 5 minutes
SLEEP_TIME=$(($RTOE * 5 * 60))
echo "NCBI BLAST search started. This may take a while (estimated time: $RTOE minutes)."
sleep $SLEEP_TIME

# Fetch the results
RESULT=$(curl -s "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=$RID")

# Output the result to a file
echo "$RESULT" > blast_report.txt

echo "BLAST search completed. Results are in blast_report.txt"
