
export MKL_NUM_THREADS=24
export NUMEXPR_NUM_THREADS=24
export OMP_NUM_THREADS=24

SCRIPTS_DIR=~/mutspec-utils/scripts

organism=human
gene=cytb


indir=data/exposure/${organism}_${gene}
echo "Input directory: $indir"

mkdir -p $indir/ms $indir/pyvolve/out

method=iqtree
replics=10
GENCODE=2

raw_tree=$indir/iqtree_anc_tree.nwk
raw_mulal=$indir/alignment_checked.fasta
spectra=$indir/ms/ms12syn.tsv

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
log_file=$indir/pyvolve/pyvolve_${method}.log

if [ ! -f ${indir}/observed_mutations_iqtree.tsv ] || 
	[ ! -f ${indir}/exp_muts_invariant.tsv ] || 
		[ ! -f $raw_mulal ] || [ ! -f $raw_tree ] ||
			[ ! -f $rates ]; then
	echo not all files are exist
	exit 1
fi
echo All input files are exits. Start computing

echo "Processing spectra calculation..."
python3 $SCRIPTS_DIR/calculate_mutspec.py -b ${indir}/observed_mutations_iqtree.tsv -e ${indir}/exp_muts_invariant.tsv -o $indir/ms \
    --exclude OUTGRP,ROOT --proba  --syn --plot -x pdf

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
echo > $log_file
cat ${mulal}.ingroup | perl -e '$p=$s="777"; while (<STDIN>) {chomp; if ($_=~/^>/) {$h=$_; if ($s!~/[^ACGT]/i) {print "$p\n$s\n"} $p=$h; $s="777"} else {$s=$_}} if ($s!~/[^ACGT]/i) {print "$p\n$s\n"}' > ${mulal}.clean
if [ `grep -c ">" ${mulal}.clean` -lt 1 ]; then
	echo -e "There are no sequences without ambigous nucleotides in the alignment.\nInterrupted" > $log_file
	exit 0
fi

python3 $SCRIPTS_DIR/pyvolve_process.py -a ${mulal}.clean -t ${tree}.ingroup -s $spectra \
	-o $indir/pyvolve/seqfile.fasta -r $replics -c $GENCODE -l 1 --write_anc --rates $rates
echo -e "Mutation samples generated\n"


for fasta_file in $indir/pyvolve/seqfile_sample-*.fasta; do
	echo "Processing $fasta_file"
	python3 $SCRIPTS_DIR/alignment2iqtree_states.py $fasta_file $fasta_file.state
	python3 $SCRIPTS_DIR/3.collect_mutations.py --tree ${tree}.ingroup --states ${fasta_file}.state \
		--gencode $GENCODE --syn --no-mutspec --outdir $indir/pyvolve/mout --force
	rm $fasta_file ${fasta_file}.state
    cat $indir/pyvolve/mout/run.log >> $log_file
	echo -e "\n\n">> $log_file
	cat $indir/pyvolve/mout/mutations.tsv > $fasta_file.mutations
done

python3 $SCRIPTS_DIR/concat_mutations.py $indir/pyvolve/seqfile_sample-*.fasta.mutations $indir/pyvolve/mutations_${method}_pyvolve.tsv
rm $indir/pyvolve/seqfile_sample-*.fasta.mutations
echo "Mutations concatenation done"

python3 $SCRIPTS_DIR/calculate_mutspec.py -b $indir/pyvolve/mutations_${method}_pyvolve.tsv -e $indir/exp_muts_invariant.tsv \
	-o $indir/pyvolve/out -l $gene --exclude OUTGRP,ROOT --syn --mnum192 0 --plot -x pdf \
    --substract192 $indir/ms/ms192syn.tsv --substract12 $indir/ms/ms12syn.tsv
echo "Mutational spectrum calculated"
