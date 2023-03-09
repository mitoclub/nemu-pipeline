THREADS = 4

if (!params.sequence){params.sequence = ""} 
if (!params.gencode){params.gencode = ""} 
if (!params.nspecies){params.nspecies = ""} 
if (!params.outgroup){params.outgroup = ""} 

params.outdir = params.sequence.replaceFirst(/\.fasta/, "")

g_2_multipleFasta_g_405 = file(params.sequence, type: 'any') 
Channel.value(params.gencode).into{g_396_gencode_g_410;g_396_gencode_g_411;g_396_gencode_g_412;g_396_gencode_g_413}
Channel.value(params.nspecies).into{g_397_mode_g_410;g_397_mode_g_411}
Channel.value(params.outgroup).set{g_398_outgroup_g_405}


process mulal_qc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /alignment_checked.fasta$/) "sequences/$filename"
	else if (filename =~ /species_mapping.txt$/) "sequences/$filename"
}

input:
 file query from g_2_multipleFasta_g_405
 val outgrp from g_398_outgroup_g_405

output:
 set val("alignment_checked"),file("alignment_checked.fasta")  into g_405_nucl_mulal_g_70, g_405_nucl_mulal_g_129, g_405_nucl_mulal_g_151, g_405_nucl_mulal_g_412, g_405_nucl_mulal_g_413
 file "species_mapping.txt" optional true  into g_405_outputFileTxt
 file "char_numbers.log"  into g_405_logFile

"""
if [ `grep -c ">" $query` -lt 6 ]; then
	echo "Number of sequences must be > 5"
	exit 0
fi

if [ `grep -E -c ">$outgrp" $query` -ne 1 ]; then
	echo "Cannot fing outgroup header in the alignment"
	exit 1
fi

grep -v  ">" $query | grep -o . | sort | uniq -c | sort -nr > char_numbers.log
if [ `head -n 5 char_numbers.log | grep -Ec "[ACGTacgt]"` -ge 3 ] && [ `grep -Ec "[EFILPQU]" char_numbers.log` -eq 0 ]; then
	echo "All right"
else
	echo "Query fasta must contain nucleotides"
	exit 1
fi

multifasta_coding.py -a $query -g $outgrp -o alignment_coded.fasta -m species_mapping.txt

/opt/dolphin/scripts/macse2.pl alignment_coded.fasta alignment_checked.fasta

"""
}


process terminal_genomes_states {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /leaves_states.state$/) "tmp/$filename"
}

input:
 set val(name), file(mulal) from g_405_nucl_mulal_g_151

output:
 set val("leaves_states"), file("leaves_states.state")  into g_151_state_g_410, g_151_state_g_411

"""
alignment2iqtree_states.py $mulal leaves_states.state
"""
}


process convert_alignment_to_phylip {

input:
 set val(name), file(mulal) from g_405_nucl_mulal_g_129

output:
 set val(name), file("${name}.phy")  into g_129_phylip_g_130, g_129_phylip_g_326, g_129_phylip_g_409, g_129_phylip_g_414

"""
java -jar /opt/readseq.jar -a -f Phylip -o ${name}.phy $mulal

"""
}

run_RAXML = params.raxml_build_tree.run_RAXML
raxml_model = params.raxml_build_tree.raxml_model
//* @style @condition:{run_RAXML="true", raxml_model}

params.RAxML_model = raxml_model


process raxml_build_tree {

input:
 set val(name), file(mulal) from g_129_phylip_g_130

output:
 set val("raxml"), file("raxml.nwk")  into g_130_tree_g_317

when:
run_RAXML == "true"

script:
"""
/opt/standard-RAxML-8.2.12/raxmlHPC-PTHREADS-SSE3 -x 987654 -p 987654 -T $THREADS -N 50 -f a -s $mulal -n Rax_tree -m $raxml_model
mv RAxML_bestTree.Rax_tree raxml.nwk
"""

}


process shrink_raxml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 set val(name), file(tree) from g_130_tree_g_317

output:
 set val("${name}_shrinked"), file("${name}_shrinked.nwk")  into g_317_tree_g_301, g_317_tree_g_109
 file "*.log"  into g_317_logFile

"""
if [ `nw_stats $tree | grep leaves | cut -f 2` -gt 8 ]; then
	run_treeshrink.py -t $tree -O treeshrink -o . -q $params.quantile -x OUTGRP
	mv treeshrink.nwk ${name}_shrinked.nwk
	mv treeshrink_summary.txt ${name}_treeshrink.log
	mv treeshrink.txt ${name}_pruned_nodes.log
else
	cat $tree > ${name}_shrinked.nwk
	echo "Shrinking are useless on such a small number of sequences" > ${name}_treeshrink.log
fi
"""
}


process extract_terminal_branch_lengths_raxml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.branches$/) "RAxML/$filename"
}

input:
 set val(name), file(tree) from g_317_tree_g_109

output:
 set val(name), file("${name}.branches")  into g_109_branches

"""
nw_distance -m p -s l -n $tree | sort -grk 2 1> ${name}.branches

if [ `grep OUTGRP ${name}.branches | cut -f 2 | python3 -c "import sys; print(float(sys.stdin.readline().strip()) > 0)"` == False ]; then
	cat "${name}.branches"
	echo "Something went wrong: outgroup is not furthest leaf in the tree"
	exit 1
fi

"""
}


process rooting_raxml_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.nwk$/) "tmp/$filename"
}

input:
 set val(name), file(tree) from g_317_tree_g_301

output:
 set val("${name}_rooted"), file("*.nwk")  into g_301_tree_g_414

"""
nw_reroot -l $tree OUTGRP 1>${name}_rooted.nwk

"""
}


process raxml_anc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /RAxML_nodeLabelledRootedTree.nwk$/) "RAxML/$filename"
	else if (filename =~ /RAxML_marginalAncestralProbabilities.state$/) "RAxML/$filename"
	else if (filename =~ /RAxML_anc_rec.log$/) "logs/$filename"
	else if (filename =~ /raw_raxml_states.txt$/) "tmp/$filename"
}

input:
 set val(namet), file(tree) from g_301_tree_g_414
 set val(name), file(mulal) from g_129_phylip_g_414

output:
 set val("raxml"),file("RAxML_nodeLabelledRootedTree.nwk")  into g_414_tree_g_411, g_414_tree_g_413
 set val("raxml"), file("RAxML_marginalAncestralProbabilities.state")  into g_414_state_g_411
 file "RAxML_marginalAncestralStates.fasta"  into g_414_multipleFasta
 file "RAxML_anc_rec.log"  into g_414_logFile
 file "raw_raxml_states.txt" optional true  into g_414_outputFileTxt

"""
nw_labels -I $tree > leaves.txt
filter_aln.py -a $mulal -l leaves.txt -o filtered_aln.phy

raxmlHPC-PTHREADS-SSE3 -T $THREADS -f A -m $params.RAxML_model -s filtered_aln.phy -t $tree -n ANCESTORS
mv RAxML_info.ANCESTORS RAxML_anc_rec.log
mv RAxML_marginalAncestralStates.ANCESTORS RAxML_marginalAncestralStates.fasta
# mv RAxML_marginalAncestralProbabilities.ANCESTORS RAxML_marginalAncestralProbabilities.txt
# mv RAxML_nodeLabelledRootedTree.ANCESTORS RAxML_nodeLabelledRootedTree.nwk

raxml_states2iqtree_states.py RAxML_marginalAncestralProbabilities.ANCESTORS RAxML_marginalAncestralProbabilities.state
rename_internal_nodes.py $tree RAxML_nodeLabelledRootedTree.ANCESTORS RAxML_nodeLabelledRootedTree.nwk


mv RAxML_marginalAncestralProbabilities.ANCESTORS raw_raxml_states.txt
"""
}


process mutations_raxml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "mutspec_tables/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
	else if (filename =~ /.*.pdf$/) "mutspec_images/$filename"
}

input:
 set val(namet), file(tree) from g_414_tree_g_411
 set val(label), file(states1) from g_414_state_g_411
 set val(names2), file(states2) from g_151_state_g_411
 val nspecies from g_397_mode_g_411
 val gencode from g_396_gencode_g_411

output:
 file "*.tsv"  into g_411_outputFileTSV
 file "*.log"  into g_411_logFile
 file "*.pdf"  into g_411_outputFilePdf
 file "ms*syn_${label}.txt" optional true  into g_411_outputFileTxt_g_413

"""
if [ $nspecies == "single" ]; then
    collect_mutations.py --tree $tree --states $states1 --states $states2 \
        --gencode $gencode --syn $params.syn4f $params.proba --no-mutspec --outdir mout
	
	mv mout/* .
    mv mutations.tsv observed_mutations_${label}.tsv
    mv run.log ${label}_mut_extraction.log

    calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_mutations.tsv -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 -l $label $params.proba --proba_min $params.proba_min --syn $params.syn4f $params.all --plot -x pdf
    calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_mutations.tsv -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 -l $label $params.proba --proba_min $params.proba_min --syn $params.syn4f $params.all --plot -x pdf --subset internal

    cp ms12syn_${label}.tsv ms12syn_${label}.txt
    cp ms192syn_${label}.tsv ms192syn_${label}.txt

elif [ $nspecies == "multiple" ]; then
    collect_mutations.py --tree $tree --states $states1 --states $states2 \
        --gencode $gencode --syn $params.syn4f $params.proba --outdir mout
        
	mv mout/* .
    mv mutations.tsv observed_mutations_${label}.tsv
    mv run.log ${label}_mut_extraction.log

    plot_spectra.py -s mutspec12.tsv  -o ${label}_ms12.pdf  -t
    plot_spectra.py -s mutspec192.tsv -o ${label}_ms192.pdf -t

else
    echo "ArgumentError: nspecies"
    exit 1
fi

"""
}


process pyvolve_raxml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "simulation/$filename"
	else if (filename =~ /.*.pdf$/) "mutspec_images/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 file spectra from g_411_outputFileTxt_g_413
 set val(label), file(tree) from g_414_tree_g_413
 set val(name), file(mulal) from g_405_nucl_mulal_g_413
 val gencode from g_396_gencode_g_413

output:
 file "*.tsv"  into g_413_outputFileTSV
 file "*.pdf"  into g_413_outputFilePdf
 file "*.log"  into g_413_logFile

when:
params.run_pyvolve == "true" && params.nspecies == "single"

script:
"""
arrray=($spectra)

spectra12=\${arrray[0]}
spectra192=\${arrray[1]}

if [-f \$spectra192 ]; then
	substractor192="--substract192 \$spectra192"
else
	substractor192=""
fi

nw_prune $tree OUTGRP | python3 /home/dolphin/dolphin/scripts/resci.py > ${tree}.ingroup
echo "Tree outgroup pruned"
awk '/^>/ {P=index(\$1, "OUTGRP")==0} {if(P) print}' $mulal > ${mulal}.ingroup
echo "Tree outgroup sequence filtered out"

nw_labels -L ${tree}.ingroup | grep -v ROOT | xargs -I{} echo -e "{}\tNode{}" > map.txt
if [ `grep -c NodeNode map.txt` -eq 0 ]; then
	nw_rename ${tree}.ingroup map.txt > ${tree}.tmp
	cat ${tree}.tmp > ${tree}.ingroup
	echo "Internal nodes renamed"
fi

#filter out sequences with ambigous nucleotides
cat ${mulal}.ingroup | perl -e '\$p=\$s="777"; while (<STDIN>) {chomp; if (\$_=~/^>/) {\$h=\$_; if (\$s!~/[^ACGT]/i) {print "\$p\n\$s\n"} \$p=\$h; \$s="777"} else {\$s=\$_}} if (\$s!~/[^ACGT]/i) {print "\$p\n\$s\n"}' > ${mulal}.clean
cat ${mulal}.clean 
if [ `grep -c ">" ${mulal}.clean` -lt 1 ]; then
	echo -e "There are no sequences without ambigous nucleotides in the alignment.\nInterrupted" > pyvolve_${label}.log
	exit 0
fi

pyvolve_process.py -a ${mulal}.clean -t ${tree}.ingroup -s \$spectra12 -o seqfile.fasta -r $params.replics -c $gencode -l $params.scale_tree --write_anc
echo "Mutation samples generated"

for fasta_file in seqfile_sample-*.fasta
do
	echo "Processing \$fasta_file"
	alignment2iqtree_states.py \$fasta_file  \${fasta_file}.state
	3.collect_mutations.py --tree ${tree}.ingroup --states  \${fasta_file}.state --gencode $gencode --syn $params.syn4f --no-mutspec --outdir mout --force
	cat mout/run.log >> pyvolve_${label}.log
	echo -e "\n\n">> pyvolve_${label}.log
	cat mout/mutations.tsv >  \${fasta_file}.mutations
done
echo "Mutations extraction done"

concat_mutations.py seqfile_sample-*.fasta.mutations mutations_${label}_pyvolve.tsv
echo "Mutations concatenation done"
cat mutations_${label}_pyvolve.tsv

calculate_mutspec.py -b mutations_${label}_pyvolve.tsv -e mout/expected_mutations.tsv -o . \
	-l ${label}_simulated --exclude OUTGRP,ROOT --syn $params.syn4f $params.all --mnum192 $params.mnum192 --plot -x pdf \
	--substract12 \$spectra12 \$substractor192
echo "Mutational spectrum calculated"

"""
}

run_IQTREE = params.iqtree_build_tree.run_IQTREE
iqtree_model = params.iqtree_build_tree.iqtree_model
quantile = params.iqtree_build_tree.quantile
//* @style @condition:{run_IQTREE="true", iqtree_model, quantile}

params.IQTREE_model = iqtree_model
params.quantile = quantile

process iqtree_build_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 set val(name), file(mulal) from g_129_phylip_g_409

output:
 set val("iqtree"), file("iqtree.nwk")  into g_409_tree_g_315
 file "*.log"  into g_409_logFile

when:
run_IQTREE == "true"

errorStrategy 'retry'
maxRetries 3

script:

"""
iqtree2 -s $mulal -m $iqtree_model -nt $THREADS --prefix phylo
mv phylo.treefile iqtree.nwk
mv phylo.iqtree iqtree_report.log
mv phylo.log iqtree.log
"""

}


process shrink_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 set val(name), file(tree) from g_409_tree_g_315

output:
 set val("${name}_shrinked"), file("${name}_shrinked.nwk")  into g_315_tree_g_302, g_315_tree_g_132
 file "*.log"  into g_315_logFile

"""
if [ `nw_stats $tree | grep leaves | cut -f 2` -gt 8 ]; then
	run_treeshrink.py -t $tree -O treeshrink -o . -q $params.quantile -x OUTGRP
	mv treeshrink.nwk ${name}_shrinked.nwk
	mv treeshrink_summary.txt ${name}_treeshrink.log
	mv treeshrink.txt ${name}_pruned_nodes.log
else
	cat $tree > ${name}_shrinked.nwk
	echo "Shrinking are useless on such a small number of sequences" > ${name}_treeshrink.log
fi
"""
}


process extract_terminal_branch_lengths_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.branches$/) "IQTREE/$filename"
}

input:
 set val(name), file(tree) from g_315_tree_g_132

output:
 set val(name), file("${name}.branches")  into g_132_branches

"""
nw_distance -m p -s l -n $tree | sort -grk 2 1> ${name}.branches

if [ `grep OUTGRP ${name}.branches | cut -f 2 | python3 -c "import sys; print(float(sys.stdin.readline().strip()) > 0)"` == False ]; then
	cat "${name}.branches"
	echo "Something went wrong: outgroup is not furthest leaf in the tree"
	exit 1
fi

"""
}


process rooting_iqtree_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.nwk$/) "tmp/$filename"
}

input:
 set val(name), file(tree) from g_315_tree_g_302

output:
 set val("${name}_rooted"), file("*.nwk")  into g_302_tree_g_326

"""
nw_reroot -l $tree OUTGRP 1>${name}_rooted.nwk

"""
}


process iqtree_anc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /iqtree_anc_tree.nwk$/) "IQTREE/$filename"
	else if (filename =~ /iqtree_anc.state$/) "IQTREE/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 set val(name), file(mulal) from g_129_phylip_g_326
 set val(namet), file(tree) from g_302_tree_g_326

output:
 set val("iqtree"), file("iqtree_anc_tree.nwk")  into g_326_tree_g_410, g_326_tree_g_412
 set val("iqtree"), file("iqtree_anc.state")  into g_326_state_g_410
 file "*.log"  into g_326_logFile

errorStrategy 'retry'
maxRetries 3

script:
"""
iqtree2 -s $mulal -m $params.IQTREE_model -asr -nt $THREADS --prefix anc
mv anc.iqtree iqtree_anc_report.log
# mv anc.state iqtree_anc.state
mv anc.log iqtree_anc.log
nw_reroot anc.treefile OUTGRP | sed 's/;/ROOT;/' > iqtree_anc_tree.nwk

iqtree_states_add_part.py anc.state iqtree_anc.state
"""

}


process mutations_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "mutspec_tables/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
	else if (filename =~ /.*.pdf$/) "mutspec_images/$filename"
}

input:
 set val(namet), file(tree) from g_326_tree_g_410
 set val(label), file(states1) from g_326_state_g_410
 set val(names2), file(states2) from g_151_state_g_410
 val nspecies from g_397_mode_g_410
 val gencode from g_396_gencode_g_410

output:
 file "*.tsv"  into g_410_outputFileTSV
 file "*.log"  into g_410_logFile
 file "*.pdf"  into g_410_outputFilePdf
 file "ms*syn_${label}.txt" optional true  into g_410_outputFileTxt_g_412

"""
if [ $nspecies == "single" ]; then
    collect_mutations.py --tree $tree --states $states1 --states $states2 \
        --gencode $gencode --syn $params.syn4f $params.proba --no-mutspec --outdir mout
	
	mv mout/* .
    mv mutations.tsv observed_mutations_${label}.tsv
    mv run.log ${label}_mut_extraction.log

    calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_mutations.tsv -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 -l $label $params.proba --proba_min $params.proba_min --syn $params.syn4f $params.all --plot -x pdf
    calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_mutations.tsv -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 -l $label $params.proba --proba_min $params.proba_min --syn $params.syn4f $params.all --plot -x pdf --subset internal

    cp ms12syn_${label}.tsv ms12syn_${label}.txt
    cp ms192syn_${label}.tsv ms192syn_${label}.txt

elif [ $nspecies == "multiple" ]; then
    collect_mutations.py --tree $tree --states $states1 --states $states2 \
        --gencode $gencode --syn $params.syn4f $params.proba --outdir mout
        
	mv mout/* .
    mv mutations.tsv observed_mutations_${label}.tsv
    mv run.log ${label}_mut_extraction.log

    plot_spectra.py -s mutspec12.tsv  -o ${label}_ms12.pdf  -t
    plot_spectra.py -s mutspec192.tsv -o ${label}_ms192.pdf -t

else
    echo "ArgumentError: nspecies"
    exit 1
fi

"""
}


process pyvolve_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "simulation/$filename"
	else if (filename =~ /.*.pdf$/) "mutspec_images/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 file spectra from g_410_outputFileTxt_g_412
 set val(label), file(tree) from g_326_tree_g_412
 set val(name), file(mulal) from g_405_nucl_mulal_g_412
 val gencode from g_396_gencode_g_412

output:
 file "*.tsv"  into g_412_outputFileTSV
 file "*.pdf"  into g_412_outputFilePdf
 file "*.log"  into g_412_logFile

when:
params.run_pyvolve == "true" && params.nspecies == "single"

script:
"""
arrray=($spectra)

spectra12=\${arrray[0]}
spectra192=\${arrray[1]}

if [-f \$spectra192 ]; then
	substractor192="--substract192 \$spectra192"
else
	substractor192=""
fi

nw_prune $tree OUTGRP | python3 /home/dolphin/dolphin/scripts/resci.py > ${tree}.ingroup
echo "Tree outgroup pruned"
awk '/^>/ {P=index(\$1, "OUTGRP")==0} {if(P) print}' $mulal > ${mulal}.ingroup
echo "Tree outgroup sequence filtered out"

nw_labels -L ${tree}.ingroup | grep -v ROOT | xargs -I{} echo -e "{}\tNode{}" > map.txt
if [ `grep -c NodeNode map.txt` -eq 0 ]; then
	nw_rename ${tree}.ingroup map.txt > ${tree}.tmp
	cat ${tree}.tmp > ${tree}.ingroup
	echo "Internal nodes renamed"
fi

#filter out sequences with ambigous nucleotides
cat ${mulal}.ingroup | perl -e '\$p=\$s="777"; while (<STDIN>) {chomp; if (\$_=~/^>/) {\$h=\$_; if (\$s!~/[^ACGT]/i) {print "\$p\n\$s\n"} \$p=\$h; \$s="777"} else {\$s=\$_}} if (\$s!~/[^ACGT]/i) {print "\$p\n\$s\n"}' > ${mulal}.clean
cat ${mulal}.clean 
if [ `grep -c ">" ${mulal}.clean` -lt 1 ]; then
	echo -e "There are no sequences without ambigous nucleotides in the alignment.\nInterrupted" > pyvolve_${label}.log
	exit 0
fi

pyvolve_process.py -a ${mulal}.clean -t ${tree}.ingroup -s \$spectra12 -o seqfile.fasta -r $params.replics -c $gencode -l $params.scale_tree --write_anc
echo "Mutation samples generated"

for fasta_file in seqfile_sample-*.fasta
do
	echo "Processing \$fasta_file"
	alignment2iqtree_states.py \$fasta_file  \${fasta_file}.state
	3.collect_mutations.py --tree ${tree}.ingroup --states  \${fasta_file}.state --gencode $gencode --syn $params.syn4f --no-mutspec --outdir mout --force
	cat mout/run.log >> pyvolve_${label}.log
	echo -e "\n\n">> pyvolve_${label}.log
	cat mout/mutations.tsv >  \${fasta_file}.mutations
done
echo "Mutations extraction done"

concat_mutations.py seqfile_sample-*.fasta.mutations mutations_${label}_pyvolve.tsv
echo "Mutations concatenation done"
cat mutations_${label}_pyvolve.tsv

calculate_mutspec.py -b mutations_${label}_pyvolve.tsv -e mout/expected_mutations.tsv -o . \
	-l ${label}_simulated --exclude OUTGRP,ROOT --syn $params.syn4f $params.all --mnum192 $params.mnum192 --plot -x pdf \
	--substract12 \$spectra12 \$substractor192
echo "Mutational spectrum calculated"

"""
}


process count_prior_mutations {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /mutnumbers.tsv$/) "logs/$filename"
}

input:
 set val(name), file(seqs) from g_405_nucl_mulal_g_70

output:
 file "mutnumbers.tsv"  into g_70_outputFileTSV

"""
/opt/dolphin/scripts/mutnumbers.pl $seqs 1>mutnumbers.tsv

"""
}


process mut_processing_params {



syn4f = params.mut_processing_params.syn4f
all = params.mut_processing_params.all
mnum192 = params.mut_processing_params.mnum192
use_probabilities = params.mut_processing_params.use_probabilities
proba_min = params.mut_processing_params.proba_min
run_pyvolve = params.mut_processing_params.run_pyvolve
replics = params.mut_processing_params.replics
scale_tree = params.mut_processing_params.scale_tree
//* @style @multicolumn:{syn4f, all, mnum192}, {use_probabilities, proba_min}, {run_pyvolve, replics, scale_tree} @condition:{run_pyvolve="true", replics, scale_tree}

println ""
println "Arguments:"
println "syn4f: ${syn4f}"
println "all: ${all}"
println "mnum192: ${mnum192}"
println "use_probabilities: ${use_probabilities}"
println "proba_min: ${proba_min}"
println "run_pyvolve: ${run_pyvolve}"
println "replics: ${replics}"
println "scale_tree: ${scale_tree}"
println ""

params.syn4f = syn4f == "true" ? "--syn4f" : ""
params.all = all == "true" ? "--all" : ""
params.mnum192 = mnum192
params.proba = use_probabilities == "true" ? "--proba" : ""
params.proba_min = proba_min
params.run_pyvolve = run_pyvolve
params.replics = replics
params.scale_tree = scale_tree

script:
"""
echo Nothing
"""


}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}