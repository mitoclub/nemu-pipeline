#!/bin/bash

# mkdir -p /home/kpotoh/nemu-pipeline/tmp/evolve_sc_raxml
cd /home/mr/nemu-pipeline/tmp/evolve

raw_tree=../evolve/iqtree_anc_tree.nwk
spectra=../evolve/ms12syn_iqtree.tsv
mulal=../evolve/alignment_checked.fasta

SCRIPTS_DIR=/home/mr/mutspec-utils/scripts

label=iqtree
replics=10
GENCODE=2
scale_tree=1

# export MKL_NUM_THREADS=1
# export NUMEXPR_NUM_THREADS=1
# export OMP_NUM_THREADS=1

tree=tree.nwk
nw_prune $raw_tree OUTGRP > ${tree}.ingroup
echo "Tree outgroup pruned"
awk '/^>/ {P=index(\$1, "OUTGRP")==0} {if(P) print}' $mulal > ${mulal}.ingroup
echo "Tree outgroup sequence filtered out"


nw_labels -L ${tree}.ingroup | grep -v ROOT | xargs -I{} echo -e "{}\tNode{}" > map.tsv
if [ `grep -c NodeNode map.tsv` -eq 0 ]; then
	nw_rename ${tree}.ingroup map.txt > ${tree}.tmp
	cat ${tree}.tmp > ${tree}.ingroup
	echo "Internal nodes renamed"
fi

#filter out sequences with ambigous nucleotides
cat ${mulal}.ingroup | perl -e '\$p=\$s=""; while (<STDIN>) {chomp; if (\$_=~/^>/) {\$h=\$_; if (\$s!~/[^ACGT]/i) {print "\$p\n\$s\n"} \$p=\$h; \$s=""} else {\$s.=\$_}} if (\$s!~/[^ACGT]/i) {print "\$p\n\$s\n"}' > ${mulal}.clean
if [ `grep -c ">" ${mulal}.clean` -lt 1 ]; then
	echo -e "There are no sequences without ambigous nucleotides in the alignment.\nInterrupted" > pyvolve_${label}.log
	exit 0
fi

python3 $SCRIPTS_DIR/pyvolve_process.py -a ${mulal}.clean -t ${tree}.ingroup -s $spectra -o seqfile.fasta -r $replics -c $GENCODE -l $scale_tree --write_anc
echo "Mutation samples generated"

echo > pyvolve_full_$label.log

for fasta_file in seqfile_sample-*.fasta; do
	echo "Processing \$fasta_file"
	python3 $SCRIPTS_DIR/alignment2iqtree_states.py $fasta_file $fasta_file.state
	python3 $SCRIPTS_DIR/3.collect_mutations.py --tree ${tree}.ingroup --states ${fasta_file}.state --gencode $GENCODE --syn --no-mutspec --outdir mout --force
	cat mout/run.log >> pyvolve_full_$label.log
	echo -e "\n\n">> pyvolve_full_$label.log
	cat mout/mutations.tsv > $fasta_file.mutations
done

python3 $SCRIPTS_DIR/concat_mutations.py seqfile_sample-*.fasta.mutations mutations_${label}_pyvolve.tsv
echo "Mutations concatenation done"

mkdir -p out
python3 $SCRIPTS_DIR/calculate_mutspec.py -b mutations_${label}_pyvolve.tsv -e mout/expected_mutations.tsv \
	-o out -l debug --exclude OUTGRP,ROOT --syn --mnum192 0 --plot -x pdf
echo "Mutational spectrum calculated"
