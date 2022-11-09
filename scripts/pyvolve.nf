# input:
#     file spectra
#     set val(label), file(tree)
#     set val(name), file(mulal)

# output:
#     file("*.tsv")
#     file("*.pdf")
#     file("*.log")

params.OUTGRP=OUTGRP
params.gencode=2
replics=10
scale_tree=1
mnum192=16
syn4f=--syn4f

nw_prune $tree $params.OUTGRP > ${tree}.ingroup
python3 pyvolve_process.py -a $mulal -t ${tree}.ingroup -s $spectra -o seqfile.fasta -r $replics -c $params.gencode -l $scale_tree

for fasta_file in seqfile_sample-*.fasta
do
	alignment2iqtree_states.py \$fasta_file  \${fasta_file}.state
	3.collect_mutations.py --tree ${tree}.ingroup --states  \${fasta_file}.state --gencode $params.gencode --syn $syn4f --no-mutspec --outdir mout --force
	cat mout/run.log >> pyvolve_${label}.log
	echo -e "\n\n">> pyvolve_${label}.log
	cat mout/mutations.tsv >  \${fasta_file}.mutations
done

concat_mutations.py seqfile_sample-*.fasta.mutations mutations_${label}_pyvolve.tsv

calculate_mutspec.py -b mutations_${label}_pyvolve.tsv -e mout/expected_mutations.tsv -o . \
	-l ${label}_pyvolve --exclude ${params.OUTGRP},ROOT --syn $syn4f --mnum192 $mnum192 --plot -x pdf
