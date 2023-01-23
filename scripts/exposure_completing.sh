
SCRIPTS_DIR=/home/mr/mutspec-utils/scripts

indir=data/exposure/human_cytb

mkdir -p $indir/ms $indir/pyvolve

echo "Processing spectra calculation..."
python3 $SCRIPTS_DIR/calculate_mutspec.py -b ${indir}/observed_mutations_iqtree.tsv -e ${indir}/exp_muts_invariant.tsv -o $indir/ms \
    --exclude OUTGRP,ROOT --proba  --syn --plot -x pdf


label=iqtree
replics=10
GENCODE=2
scale_tree=1

raw_tree=$indir/iqtree_anc_tree.nwk
spectra=$indir/ms/ms12syn.tsv
raw_mulal=$indir/alignment_checked.fasta

tree=$indir/pyvolve/tree.nwk
mulal=$indir/pyvolve/mulal.fasta
log_file=$indir/pyvolve/pyvolve_${label}.log

# export MKL_NUM_THREADS=1
# export NUMEXPR_NUM_THREADS=1
# export OMP_NUM_THREADS=1


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

python3 $SCRIPTS_DIR/pyvolve_process.py -a ${mulal}.clean -t ${tree}.ingroup -s $spectra -o $indir/pyvolve/seqfile.fasta -r $replics -c $GENCODE -l $scale_tree --write_anc
echo "Mutation samples generated"


for fasta_file in $indir/pyvolve/seqfile_sample-*.fasta; do
	echo "Processing $fasta_file"
	python3 $SCRIPTS_DIR/alignment2iqtree_states.py $fasta_file $fasta_file.state
	python3 $SCRIPTS_DIR/3.collect_mutations.py --tree ${tree}.ingroup --states ${fasta_file}.state --gencode $GENCODE --syn --no-mutspec --outdir $indir/pyvolve/mout --force
	rm ${fasta_file}.state
    cat $indir/pyvolve/mout/run.log >> $log_file
	echo -e "\n\n">> $log_file
	cat $indir/pyvolve/mout/mutations.tsv > $fasta_file.mutations
done

python3 $SCRIPTS_DIR/concat_mutations.py $indir/pyvolve/seqfile_sample-*.fasta.mutations $indir/pyvolve/mutations_${label}_pyvolve.tsv
echo "Mutations concatenation done"

mkdir -p $indir/pyvolve/out
python3 $SCRIPTS_DIR/calculate_mutspec.py -b $indir/pyvolve/mutations_${label}_pyvolve.tsv -e $indir/exp_muts_invariant.tsv \
	-o $indir/pyvolve/out -l debug --exclude OUTGRP,ROOT --syn --mnum192 0 --plot -x pdf
echo "Mutational spectrum calculated"
