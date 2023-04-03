THREADS=128
export MKL_NUM_THREADS=$THREADS
export NUMEXPR_NUM_THREADS=$THREADS
export OMP_NUM_THREADS=$THREADS

# organism=mus
# gene=nd1
for organism in mus human; do
for gene in nd1 cytb; do


################################################################################
################################ PREPARATION ###################################
################################################################################

indir=data/exposure/${organism}_${gene}
echo "Input directory: $indir"

mkdir -p $indir/ms $indir/pyvolve/out

method=iqtree
replics=100
GENCODE=2

raw_tree=$indir/iqtree_anc_tree.nwk
raw_mulal=$indir/alignment_checked.fasta

if [ $gene == cytb ]; then
	rates=$indir/CYTB.rate
elif [ $gene == nd1 ]; then
	rates=$indir/ND1.rate
else
	echo NotImplementedError
	exit 1
fi
# echo "rates: $rates" 

tree=$indir/pyvolve/tree.nwk
mulal=$indir/pyvolve/mulal.fasta

if [ ! -f ${indir}/observed_mutations_iqtree.tsv ] || 
	[ ! -f ${indir}/exp_muts_invariant.tsv ] || 
		[ ! -f $raw_mulal ] || [ ! -f $raw_tree ] ||
			[ ! -f $rates ]; then
	echo not all files are exist
	exit 1
fi
echo All input files are exits. Start computing

################################################################################
############################### START COMPUTING ################################
################################################################################

echo "Processing spectra calculation..."
# calculate_mutspec.py -b ${indir}/observed_mutations_iqtree.tsv -e ${indir}/exp_muts_invariant.tsv -o $indir/ms \
#     --exclude OUTGRP,ROOT --proba  --syn --plot -x pdf -l ${organism}_${gene} --mnum192 0
calculate_mutspec.py -b ${indir}/observed_mutations_iqtree.tsv -e ${indir}/exp_muts_invariant.tsv -o $indir/ms \
    --exclude OUTGRP,ROOT --proba  --syn --plot -x pdf -l ${organism}_${gene} --subset internal --mnum192 0

# tree processing
nw_prune $raw_tree OUTGRP | python3 -c "import sys,re; print(re.sub('\d+\.\d+e-\d+', lambda m: '{:.10f}'.format(float(m.group())), sys.stdin.read().strip()))" > ${tree}.ingroup
echo "Tree outgroup pruned"
## replace digit-named nodes from RAxML
nw_labels -L ${tree}.ingroup | grep -v ROOT | xargs -I{} echo -e "{}\tNode{}" > $indir/pyvolve/map.txt
if [ `grep -c NodeNode $indir/pyvolve/map.txt` -eq 0 ]; then
	nw_rename ${tree}.ingroup $indir/pyvolve/map.txt > ${tree}.tmp
	cat ${tree}.tmp > ${tree}.ingroup
	echo "Internal nodes renamed"
fi

# alignment processing
awk '/^>/ {P=index($1, "OUTGRP")==0} {if(P) print}' $raw_mulal > ${mulal}.ingroup
echo "Tree outgroup sequence filtered out"
## filter out sequences with ambigous nucleotides
cat ${mulal}.ingroup | perl -e '$p=$s="777"; while (<STDIN>) {chomp; if ($_=~/^>/) {$h=$_; if ($s!~/[^ACGT]/i) {print "$p\n$s\n"} $p=$h; $s="777"} else {$s=$_}} if ($s!~/[^ACGT]/i) {print "$p\n$s\n"}' > ${mulal}.clean
if [ `grep -c ">" ${mulal}.clean` -lt 1 ]; then
	echo -e "There are no sequences without ambigous nucleotides in the alignment.\nInterrupted"
	exit 0
fi

pyvolve_process.py -a ${mulal}.clean -t ${tree}.ingroup \
	-s $indir/ms/ms12syn_internal_${organism}_${gene}.tsv \
	-o $indir/pyvolve/seqfile.fasta -r $replics -c $GENCODE \
	-l 1 --write_anc --rates $rates
echo -e "Mutation samples generated\n"

parallel --jobs $THREADS alignment2iqtree_states.py {} {}.state ::: $indir/pyvolve/seqfile_sample-*.fasta

parallel --jobs $THREADS collect_mutations.py --tree ${tree}.ingroup --states {} --gencode $GENCODE --syn --no-mutspec \
			 --outdir $indir/pyvolve/mout/{#} --force --quiet ::: $indir/pyvolve/seqfile_sample-*.fasta.state

rm $indir/pyvolve/seqfile_sample-*.fasta $indir/pyvolve/seqfile_sample-*.fasta.state

concat_mutations.py $indir/pyvolve/mout/*/mutations.tsv $indir/pyvolve/mutations_${method}_pyvolve.tsv
echo "Mutations concatenation done"


# calculate_mutspec.py -b $indir/pyvolve/mutations_${method}_pyvolve.tsv -e $indir/exp_muts_invariant.tsv \
# 	-o $indir/pyvolve/out -l ${organism}_${gene}_simulated --exclude OUTGRP,ROOT --syn --mnum192 0 --plot -x pdf \
#     --substract192 $indir/ms/ms192syn_${organism}_${gene}.tsv --substract12 $indir/ms/ms12syn_${organism}_${gene}.tsv
calculate_mutspec.py -b $indir/pyvolve/mutations_${method}_pyvolve.tsv -e $indir/exp_muts_invariant.tsv \
	-o $indir/pyvolve/out -l ${organism}_${gene}_internal_simulated --exclude OUTGRP,ROOT --syn --mnum192 0 --plot -x pdf \
    --substract192 $indir/ms/ms192syn_internal_${organism}_${gene}.tsv --substract12 $indir/ms/ms12syn_internal_${organism}_${gene}.tsv
echo "Mutational spectrum calculated"


done
done