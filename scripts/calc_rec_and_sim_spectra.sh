# part of run_simulation.sh script

for organism in mouse human; do
for gene in nd1 cytb cox1; do


################################################################################
################################ PREPARATION ###################################
################################################################################

indir=data/selection_search/${organism}_${gene}
echo -e "\nInput directory: $indir"

mkdir -p $indir/ms $indir/pyvolve/out

method=iqtree
GENCODE=2

raw_tree=$indir/IQTREE/iqtree_anc_tree.nwk
raw_mulal=$indir/sequences/alignment_checked.fasta
obs_mutations=$indir/tables/observed_mutations_iqtree.tsv
exp_mutations=$indir/exp_muts_invariant.tsv  # TODO create
rates=data/selection_search/rates/${organism}_${gene}.rate

if [ ! -f $obs_mutations ] || [ ! -f $exp_mutations ] || 
	   [ ! -f $raw_mulal ] || [ ! -f $raw_tree ] || [ ! -f $rates ]; then
	echo not all files are exist
	exit 1
fi
echo All input files are exits. Start computing

################################################################################
############################### START COMPUTING ################################
################################################################################

tree=$indir/pyvolve/tree.nwk
mulal=$indir/pyvolve/mulal.fasta

echo "Processing spectra calculation..."
calculate_mutspec.py -b $obs_mutations -e $exp_mutations -o $indir/ms \
    --exclude OUTGRP,ROOT --proba --syn --plot -x pdf -l ${organism}_${gene} --mnum192 0
calculate_mutspec.py -b $obs_mutations -e $exp_mutations -o $indir/ms \
    --exclude OUTGRP,ROOT --proba --syn --plot -x pdf -l ${organism}_${gene} --subset internal --mnum192 0
echo "Reconstructed mutational spectrum calculated"


# ALL BRANCHES
calculate_mutspec.py -b $indir/pyvolve/mutations_${method}_pyvolve.tsv -e $exp_mutations \
	-o $indir/pyvolve/out --exclude OUTGRP,ROOT --syn --mnum192 0 --plot -x pdf \
	-l ${organism}_${gene}_simulated \
	--substract12 $indir/ms/ms12syn_${organism}_${gene}.tsv \
    --substract192 $indir/ms/ms192syn_${organism}_${gene}.tsv \

# INTERNAL BRANCHES ONLY
calculate_mutspec.py -b $indir/pyvolve/mutations_${method}_pyvolve.tsv -e $exp_mutations \
	-o $indir/pyvolve/out  --exclude OUTGRP,ROOT --syn --mnum192 0 --plot -x pdf \
	-l ${organism}_${gene}_internal_simulated --subset internal \
	--substract12 $indir/ms/ms12syn_internal_${organism}_${gene}.tsv \
    --substract192 $indir/ms/ms192syn_internal_${organism}_${gene}.tsv \

echo "Simulated mutational spectrum calculated"

done
done
