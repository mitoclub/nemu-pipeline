THREADS=64
export MKL_NUM_THREADS=$THREADS
export NUMEXPR_NUM_THREADS=$THREADS
export OMP_NUM_THREADS=$THREADS


organism=human
gene=$1
# gene=nd1
# for organism in mouse human; do
# for gene in nd1 cytb cox1; do


################################################################################
################################ PREPARATION ###################################
################################################################################

indir=data/selection_search/${organism}_${gene}
echo -e "\nInput directory: $indir"

mkdir -p $indir/ms $indir/pyvolve/out

method=iqtree
replics=100
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
    --exclude OUTGRP,ROOT --proba --syn --plot -x pdf -l ${organism}_${gene} --subset internal --mnum192 0

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


# calculate_mutspec.py -b $indir/pyvolve/mutations_${method}_pyvolve.tsv -e $exp_mutations \
# 	-o $indir/pyvolve/out -l ${organism}_${gene}_simulated --exclude OUTGRP,ROOT --syn --mnum192 0 --plot -x pdf \
#     --substract192 $indir/ms/ms192syn_${organism}_${gene}.tsv --substract12 $indir/ms/ms12syn_${organism}_${gene}.tsv
calculate_mutspec.py -b $indir/pyvolve/mutations_${method}_pyvolve.tsv -e $exp_mutations \
	-o $indir/pyvolve/out -l ${organism}_${gene}_internal_simulated --exclude OUTGRP,ROOT --syn --mnum192 0 --plot -x pdf \
    --substract192 $indir/ms/ms192syn_internal_${organism}_${gene}.tsv --substract12 $indir/ms/ms12syn_internal_${organism}_${gene}.tsv
echo "Mutational spectrum calculated"

# done
# done
