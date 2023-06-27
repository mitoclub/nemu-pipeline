// params.outdir = params.sequence.replaceFirst(/\.fasta/, "")

if (!params.sequence){params.sequence = ""} 
if (!params.gencode){params.gencode = ""} 
if (!params.nspecies){params.nspecies = ""} 
if (!params.outgroup){params.outgroup = ""} 
if (!params.aligned){params.aligned = ""} 
if (!params.treefile){params.treefile = ""} 

g_2_multipleFasta_g_428 = file(params.sequence, type: 'any') 
precalculated_tree = file(params.treefile, type: 'any') 
Channel.value(params.gencode).into{g_396_gencode_g_410;g_396_gencode_g_411;g_396_gencode_g_422;g_396_gencode_g_423;g_396_gencode_g_433}
Channel.value(params.nspecies).into{g_397_mode_g_410;g_397_mode_g_411;g_397_mode_g_422;g_397_mode_g_423}
Channel.value(params.outgroup).set{g_398_outgroup_g_428}
Channel.value(params.aligned).set{g_431_type_g_433}


process fasta_qc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /species_mapping.txt$/) "sequences/$filename"
}

input:
 file query from g_2_multipleFasta_g_428
 val outgrp from g_398_outgroup_g_428

output:
 file "sequences.fasta"  into g_428_multipleFasta_g_433
 file "species_mapping.txt" optional true  into g_428_outputFileTxt
 file "char_numbers.log"  into g_428_logFile

"""
if [ `grep -c ">" $query` -lt 10 ]; then
	echo "Number of sequences must be >= 10"
	exit 1
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

multifasta_coding.py -a $query -g $outgrp -o sequences.fasta -m species_mapping.txt

"""
}


process macse {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /alignment_checked.fasta$/) "sequences/$filename"
}

input:
 val aligned from g_431_type_g_433
 file seqs from g_428_multipleFasta_g_433
 val gencode from g_396_gencode_g_433

output:
 file "alignment_checked.fasta"  into g_433_multipleFasta_g_420, g_433_multipleFasta_g_421, g_433_multipleFasta_g_423, g_433_multipleFasta_g_422, g_433_multipleFasta_g_424
 file "aln_aa.fasta" optional true  into g_433_multipleFasta

"""
if [ $aligned == "true" ]; then
	echo "No need to align!"
	cp $seqs aln.fasta
elif [ $aligned == "false" ]; then
	java -jar /opt/macse_v2.05.jar -prog alignSequences -gc_def $gencode -out_AA aln_aa.fasta -out_NT aln.fasta -seq $seqs
else
	echo "Inappropriate values for 'aligned'"
	exit 1
fi

echo "Do quality control"
/opt/dolphin/scripts/macse2.pl aln.fasta alignment_checked.fasta

"""
}


process convert_alignment_to_phylip {

input:
 file aln from g_433_multipleFasta_g_424

output:
 set val("aln"), file("aln.phy")  into g_424_phylip_g_130, g_424_phylip_g_326, g_424_phylip_g_409, g_424_phylip_g_430

"""
java -jar /opt/readseq.jar -a -f Phylip -o aln.phy $aln

"""
}

run_iqtree = params.phylo.run_iqtree

process iqtree_build_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.log$/) "logs/$filename"
	else if (filename =~ /.*.nwk$/) "tmp/$filename"
}

input:
 set val(name), file(mulal) from g_424_phylip_g_409
 file prectree from precalculated_tree
 file labels_mapping from g_428_outputFileTxt

output:
 set val("iqtree"), file("iqtree.nwk")  into g_409_tree_g_315
 file "*.log" optional true  into g_409_logFile

when:
run_iqtree == "true"

errorStrategy 'retry'
maxRetries 3

script:

"""
if [ -f $prectree ]; then
	if [ ! -f $labels_mapping ]; then 
		echo "Something went wrong. Expected that tree must be relabeled \
		according to relabeled alignment, but file with labels map doesn't exist"
		exit 1
	fi

	awk '{print \$2 "\t" \$1}' species_mapping.txt > species2label.txt
	nw_rename $prectree species2label.txt > iqtree.nwk

else
	iqtree2 -s $mulal -m $params.phylo.iqtree_model -nt $THREADS --prefix phylo
	mv phylo.treefile iqtree.nwk
	mv phylo.iqtree iqtree_report.log
	mv phylo.log iqtree.log
fi
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
if [ $params.phylo.run_shrinking == true ] && [ `nw_stats $tree | grep leaves | cut -f 2` -gt 8 ]; then
	run_treeshrink.py -t $tree -O treeshrink -o . -q $params.phylo.quantile -x OUTGRP
	mv treeshrink.nwk ${name}_shrinked.nwk
	mv treeshrink_summary.txt ${name}_treeshrink.log
	mv treeshrink.txt ${name}_pruned_nodes.log
else
	cat $tree > ${name}_shrinked.nwk
	if [ $params.phylo.run_shrinking == true ]; then
		echo "Shrinking are useless on such a small number of sequences" > ${name}_treeshrink.log
	fi
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
nw_distance -m p -s f -n $tree | sort -grk 2 1> ${name}.branches

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


process terminal_genomes_states {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /leaves_states.state$/) "tmp/$filename"
}

input:
 file aln from g_433_multipleFasta_g_421

output:
 set val("leaves_states"), file("leaves_states.state")  into g_421_state_g_410, g_421_state_g_411

"""
alignment2iqtree_states.py $aln leaves_states.state
"""
}


rate_arg = params.siterates.exclude_cons_sites == "true" ? "--rate" : ""

process iqtree_anc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /iqtree_anc_tree.nwk$/) "IQTREE/$filename"
	else if (filename =~ /iqtree_anc.state$/) "IQTREE/$filename"
	else if (filename =~ /.*.rate$/) "IQTREE/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 set val(name), file(mulal) from g_424_phylip_g_326
 set val(namet), file(tree) from g_302_tree_g_326

output:
 set val("iqtree"), file("iqtree_anc_tree.nwk")  into g_326_tree_g_410, g_326_tree_g_422
 set val("iqtree"), file("iqtree_anc.state")  into g_326_state_g_410
 file "*.rate"  optional true into g_326_ratefile
 file "*.log"  into g_326_logFile

errorStrategy 'retry'
maxRetries 3

script:
"""
nw_labels -I $tree > leaves.txt
filter_aln.py -a $mulal -l leaves.txt -o filtered_aln.phy

iqtree2 -te $tree -s filtered_aln.phy -m $params.phylo.iqtree_anc_model -asr -nt $THREADS --prefix anc $rate_arg
mv anc.iqtree iqtree_anc_report.log
mv anc.log iqtree_anc.log
nw_reroot anc.treefile OUTGRP | sed 's/;/ROOT;/' > iqtree_anc_tree.nwk

iqtree_states_add_part.py anc.state iqtree_anc.state
"""

}


save_exp_muts = params.mut_processing_params.save_exp_mutations == "true" ? "--save-exp-muts" : ""

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
 file rates from g_326_ratefile
 set val(names2), file(states2) from g_421_state_g_410
 val nspecies from g_397_mode_g_410
 val gencode from g_396_gencode_g_410

output:
 file "*.tsv"  into g_410_outputFileTSV
 file "*.log"  into g_410_logFile
 file "*.pdf"  into g_410_outputFilePdf
 file "ms*syn_${label}.txt" optional true  into g_410_outputFileTxt_g_422

"""
if [-f $rates ]; then
	rates_arg="--rates $rates"
	echo "TODO add external rates file support here (if exist, check that it has appropriate format and replace generated $rates file)"
else
	rates_arg=""
fi


if [ $nspecies == "single" ]; then
    collect_mutations.py --tree $tree --states $states1 --states $states2 \
        --gencode $gencode --syn $params.syn4f $params.proba_arg --no-mutspec \
		\${rates_arg} --cat-cutoff $params.siterates.cons_cat_cutoff \
		--outdir mout $save_exp_muts
	
	mv mout/* .
    mv mutations.tsv observed_mutations_${label}.tsv
    mv run.log ${label}_mut_extraction.log

    calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_freqs.tsv -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 -l $label $params.proba_arg --proba_min $params.proba_cutoff --syn $params.syn4f $params.all --plot -x pdf
    calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_freqs.tsv -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 -l $label $params.proba_arg --proba_min $params.proba_cutoff --syn $params.syn4f $params.all --plot -x pdf --subset internal

    cp ms12syn_${label}.tsv ms12syn_${label}.txt
    if [-f ms192syn_${label}.tsv ]; then
	    cp ms192syn_${label}.tsv ms192syn_${label}.txt
	fi
	
elif [ $nspecies == "multiple" ]; then
    collect_mutations.py --tree $tree --states $states1 --states $states2 \
        --gencode $gencode --syn $params.syn4f $params.proba_arg --rates $rates \
		\${rates_arg} --cat-cutoff $params.siterates.cons_cat_cutoff \
		--outdir mout
        
	mv mout/* .
    mv mutations.tsv observed_mutations_${label}.tsv
    mv run.log ${label}_mut_extraction.log

    plot_spectra.py -s mutspec12.tsv  -o ${label}_ms12.pdf
    plot_spectra.py -s mutspec192.tsv -o ${label}_ms192.pdf

else
    echo "ArgumentError: nspecies"
    exit 1
fi

"""
}


process pyvolve_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "mutspec_tables/$filename"
	else if (filename =~ /.*.pdf$/) "mutspec_images/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 file spectra from g_410_outputFileTxt_g_422
 set val(label), file(tree) from g_326_tree_g_422
 file mulal from g_433_multipleFasta_g_422
 val gencode from g_396_gencode_g_422
 val nspecies from g_397_mode_g_422

output:
 file "*.tsv"  into g_422_outputFileTSV
 file "*.pdf"  into g_422_outputFilePdf
 file "*.log"  into g_422_logFile

when:
params.run_simulation == "true" && nspecies == "single"

script:
"""
arrray=($spectra)

spectra12=\${arrray[0]}
spectra192=\${arrray[1]}

if [-f \$spectra192 ]; then
	substractor192="--substract192 \${spectra192}"
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
	collect_mutations.py --tree ${tree}.ingroup --states  \${fasta_file}.state --gencode $gencode --syn $params.syn4f --no-mutspec --outdir mout --force
	cat mout/run.log >> pyvolve_${label}.log
	echo -e "\n\n">> pyvolve_${label}.log
	cat mout/mutations.tsv >  \${fasta_file}.mutations
done
echo "Mutations extraction done"

concat_mutations.py seqfile_sample-*.fasta.mutations mutations_${label}_pyvolve.tsv
echo "Mutations concatenation done"

calculate_mutspec.py -b mutations_${label}_pyvolve.tsv -e mout/expected_freqs.tsv -o . \
	-l ${label}_simulated --exclude OUTGRP,ROOT --syn $params.syn4f $params.all --mnum192 $params.mnum192 --plot -x pdf \
	--substract12 \$spectra12 \$substractor192
echo "Mutational spectrum calculated"

"""
}



process mut_processing_params {



syn4f = params.mut_processing_params.syn4f
all = params.mut_processing_params.all
mnum192 = params.mut_processing_params.mnum192
use_probabilities = params.mut_processing_params.use_probabilities
proba_cutoff = params.mut_processing_params.proba_cutoff
run_simulation = params.mut_processing_params.run_simulation
replics = params.mut_processing_params.replics
scale_tree = params.mut_processing_params.scale_tree
//* @style @multicolumn:{syn4f, all, mnum192}, {use_probabilities, proba_cutoff}, {run_simulation, replics, scale_tree} @condition:{run_simulation="true", replics, scale_tree}

println ""
println "Arguments:"
println "all: ${all}"
println "syn: true"
println "syn4f: ${syn4f}"
println "Minimal number of mutations to save 192-component spectrum (mnum192): ${mnum192}"
println "Use probabilities: ${use_probabilities}"
println "Mutation probability cutoff: ${proba_cutoff}"
println "Run simulation: ${run_simulation}"
println "Replics in simulation: ${replics}"
println "Tree scaling coefficient: ${scale_tree}"
println "Threads: ${THREADS}"
println "Run treeshrink: ${params.phylo.run_shrinking}"
println "Shrinking Quantile: ${params.phylo.quantile}"
println "Exclude conservstive sites: ${params.siterates.exclude_cons_sites}"
println ""

params.syn4f = syn4f == "true" ? "--syn4f" : ""
params.all = all == "true" ? "--all" : ""
params.mnum192 = mnum192
params.proba_arg = use_probabilities == "true" ? "--proba" : ""
params.proba_cutoff = proba_cutoff
params.run_simulation = run_simulation
params.replics = replics
params.scale_tree = scale_tree

script:
"""
echo Nothing
"""
}



// run_RAXML = params.raxml_build_tree.run_RAXML
// raxml_model = params.raxml_build_tree.raxml_model
// //* @style @condition:{run_RAXML="true", raxml_model}

// params.RAxML_model = raxml_model


// process raxml_build_tree {

// input:
//  set val(name), file(mulal) from g_424_phylip_g_130

// output:
//  set val("raxml"), file("raxml.nwk")  into g_130_tree_g_317

// when:
// run_RAXML == "true"

// script:
// """
// /opt/standard-RAxML-8.2.12/raxmlHPC-PTHREADS-SSE3 -x 987654 -p 987654 -T $THREADS -N 50 -f a -s $mulal -n Rax_tree -m $raxml_model
// mv RAxML_bestTree.Rax_tree raxml.nwk
// """

// }


// process shrink_raxml {

// publishDir params.outdir, overwrite: true, mode: 'copy',
// 	saveAs: {filename ->
// 	if (filename =~ /.*.log$/) "logs/$filename"
// }

// input:
//  set val(name), file(tree) from g_130_tree_g_317

// output:
//  set val("${name}_shrinked"), file("${name}_shrinked.nwk")  into g_317_tree_g_301, g_317_tree_g_109
//  file "*.log"  into g_317_logFile

// """
// if [ `nw_stats $tree | grep leaves | cut -f 2` -gt 8 ]; then
// 	run_treeshrink.py -t $tree -O treeshrink -o . -q $params.phylo.quantile -x OUTGRP
// 	mv treeshrink.nwk ${name}_shrinked.nwk
// 	mv treeshrink_summary.txt ${name}_treeshrink.log
// 	mv treeshrink.txt ${name}_pruned_nodes.log
// else
// 	cat $tree > ${name}_shrinked.nwk
// 	echo "Shrinking are useless on such a small number of sequences" > ${name}_treeshrink.log
// fi
// """
// }


// process extract_terminal_branch_lengths_raxml {

// publishDir params.outdir, overwrite: true, mode: 'copy',
// 	saveAs: {filename ->
// 	if (filename =~ /${name}.branches$/) "RAxML/$filename"
// }

// input:
//  set val(name), file(tree) from g_317_tree_g_109

// output:
//  set val(name), file("${name}.branches")  into g_109_branches

// """
// nw_distance -m p -s l -n $tree | sort -grk 2 1> ${name}.branches

// if [ `grep OUTGRP ${name}.branches | cut -f 2 | python3 -c "import sys; print(float(sys.stdin.readline().strip()) > 0)"` == False ]; then
// 	cat "${name}.branches"
// 	echo "Something went wrong: outgroup is not furthest leaf in the tree"
// 	exit 1
// fi

// """
// }


// process rooting_raxml_tree {

// publishDir params.outdir, overwrite: true, mode: 'copy',
// 	saveAs: {filename ->
// 	if (filename =~ /.*.nwk$/) "tmp/$filename"
// }

// input:
//  set val(name), file(tree) from g_317_tree_g_301

// output:
//  set val("${name}_rooted"), file("*.nwk")  into g_301_tree_g_430

// """
// nw_reroot -l $tree OUTGRP 1>${name}_rooted.nwk

// """
// }


// process raxml_anc {

// publishDir params.outdir, overwrite: true, mode: 'copy',
// 	saveAs: {filename ->
// 	if (filename =~ /RAxML_nodeLabelledRootedTree.nwk$/) "RAxML/$filename"
// 	else if (filename =~ /RAxML_marginalAncestralProbabilities.state$/) "RAxML/$filename"
// 	else if (filename =~ /RAxML_anc_rec.log$/) "logs/$filename"
// 	else if (filename =~ /RAxML_marginalAncestralProbabilities.ANCESTORS.txt$/) "tmp/$filename"
// }

// input:
//  set val(namet), file(tree) from g_301_tree_g_430
//  set val(name), file(mulal) from g_424_phylip_g_430

// output:
//  set val("raxml"),file("RAxML_nodeLabelledRootedTree.nwk")  into g_430_tree_g_411, g_430_tree_g_423
//  set val("raxml"), file("RAxML_marginalAncestralProbabilities.state")  into g_430_state_g_411
//  file "RAxML_marginalAncestralStates.fasta"  into g_430_multipleFasta
//  file "RAxML_anc_rec.log"  into g_430_logFile
//  file "RAxML_marginalAncestralProbabilities.ANCESTORS.txt" optional true  into g_430_outputFileTxt

// """
// nw_labels -I $tree > leaves.txt
// filter_aln.py -a $mulal -l leaves.txt -o filtered_aln.phy

// raxmlHPC-PTHREADS-SSE3 -T $THREADS -f A -m $params.RAxML_model -s filtered_aln.phy -t $tree -n ANCESTORS
// mv RAxML_info.ANCESTORS RAxML_anc_rec.log
// mv RAxML_marginalAncestralStates.ANCESTORS RAxML_marginalAncestralStates.fasta
// # mv RAxML_marginalAncestralProbabilities.ANCESTORS RAxML_marginalAncestralProbabilities.txt
// # mv RAxML_nodeLabelledRootedTree.ANCESTORS RAxML_nodeLabelledRootedTree.nwk

// raxml_states2iqtree_states.py RAxML_marginalAncestralProbabilities.ANCESTORS RAxML_marginalAncestralProbabilities.state
// rename_internal_nodes.py $tree RAxML_nodeLabelledRootedTree.ANCESTORS RAxML_nodeLabelledRootedTree.nwk


// mv RAxML_marginalAncestralProbabilities.ANCESTORS RAxML_marginalAncestralProbabilities.ANCESTORS.txt
// """
// }


// process mutations_raxml {

// publishDir params.outdir, overwrite: true, mode: 'copy',
// 	saveAs: {filename ->
// 	if (filename =~ /.*.tsv$/) "mutspec_tables/$filename"
// 	else if (filename =~ /.*.log$/) "logs/$filename"
// 	else if (filename =~ /.*.pdf$/) "mutspec_images/$filename"
// }

// input:
//  set val(namet), file(tree) from g_430_tree_g_411
//  set val(label), file(states1) from g_430_state_g_411
//  set val(names2), file(states2) from g_421_state_g_411
//  val nspecies from g_397_mode_g_411
//  val gencode from g_396_gencode_g_411

// output:
//  file "*.tsv"  into g_411_outputFileTSV
//  file "*.log"  into g_411_logFile
//  file "*.pdf"  into g_411_outputFilePdf
//  file "ms*syn_${label}.txt" optional true  into g_411_outputFileTxt_g_423

// """
// if [ $nspecies == "single" ]; then
//     collect_mutations.py --tree $tree --states $states1 --states $states2 \
//         --gencode $gencode --syn $params.syn4f $params.proba_arg --no-mutspec --outdir mout
	
// 	mv mout/* .
//     mv mutations.tsv observed_mutations_${label}.tsv
//     mv run.log ${label}_mut_extraction.log

//     calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_freqs.tsv -o . \
//         --exclude OUTGRP,ROOT --mnum192 $params.mnum192 -l $label $params.proba_arg --proba_min $params.proba_cutoff --syn $params.syn4f $params.all --plot -x pdf
//     calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_freqs.tsv -o . \
//         --exclude OUTGRP,ROOT --mnum192 $params.mnum192 -l $label $params.proba_arg --proba_min $params.proba_cutoff --syn $params.syn4f $params.all --plot -x pdf --subset internal

//     cp ms12syn_${label}.tsv ms12syn_${label}.txt
//     if [-f ms192syn_${label}.tsv ]; then
// 	    cp ms192syn_${label}.tsv ms192syn_${label}.txt
// 	fi
	
// elif [ $nspecies == "multiple" ]; then
//     collect_mutations.py --tree $tree --states $states1 --states $states2 \
//         --gencode $gencode --syn $params.syn4f $params.proba_arg --outdir mout
        
// 	mv mout/* .
//     mv mutations.tsv observed_mutations_${label}.tsv
//     mv run.log ${label}_mut_extraction.log

//     plot_spectra.py -s mutspec12.tsv  -o ${label}_ms12.pdf
//     plot_spectra.py -s mutspec192.tsv -o ${label}_ms192.pdf

// else
//     echo "ArgumentError: nspecies"
//     exit 1
// fi

// """
// }


// process pyvolve_raxml {

// publishDir params.outdir, overwrite: true, mode: 'copy',
// 	saveAs: {filename ->
// 	if (filename =~ /.*.tsv$/) "mutspec_tables/$filename"
// 	else if (filename =~ /.*.pdf$/) "mutspec_images/$filename"
// 	else if (filename =~ /.*.log$/) "logs/$filename"
// }

// input:
//  file spectra from g_411_outputFileTxt_g_423
//  set val(label), file(tree) from g_430_tree_g_423
//  file mulal from g_433_multipleFasta_g_423
//  val gencode from g_396_gencode_g_423
//  val nspecies from g_397_mode_g_423

// output:
//  file "*.tsv"  into g_423_outputFileTSV
//  file "*.pdf"  into g_423_outputFilePdf
//  file "*.log"  into g_423_logFile

// when:
// params.run_simulation == "true" && nspecies == "single"

// script:
// """
// arrray=($spectra)

// spectra12=\${arrray[0]}
// spectra192=\${arrray[1]}

// if [-f \$spectra192 ]; then
// 	substractor192="--substract192 \${spectra192}"
// else
// 	substractor192=""
// fi

// nw_prune $tree OUTGRP | python3 /home/dolphin/dolphin/scripts/resci.py > ${tree}.ingroup
// echo "Tree outgroup pruned"
// awk '/^>/ {P=index(\$1, "OUTGRP")==0} {if(P) print}' $mulal > ${mulal}.ingroup
// echo "Tree outgroup sequence filtered out"

// nw_labels -L ${tree}.ingroup | grep -v ROOT | xargs -I{} echo -e "{}\tNode{}" > map.txt
// if [ `grep -c NodeNode map.txt` -eq 0 ]; then
// 	nw_rename ${tree}.ingroup map.txt > ${tree}.tmp
// 	cat ${tree}.tmp > ${tree}.ingroup
// 	echo "Internal nodes renamed"
// fi

// #filter out sequences with ambigous nucleotides
// cat ${mulal}.ingroup | perl -e '\$p=\$s="777"; while (<STDIN>) {chomp; if (\$_=~/^>/) {\$h=\$_; if (\$s!~/[^ACGT]/i) {print "\$p\n\$s\n"} \$p=\$h; \$s="777"} else {\$s=\$_}} if (\$s!~/[^ACGT]/i) {print "\$p\n\$s\n"}' > ${mulal}.clean
// if [ `grep -c ">" ${mulal}.clean` -lt 1 ]; then
// 	echo -e "There are no sequences without ambigous nucleotides in the alignment.\nInterrupted" > pyvolve_${label}.log
// 	exit 0
// fi

// pyvolve_process.py -a ${mulal}.clean -t ${tree}.ingroup -s \$spectra12 -o seqfile.fasta -r $params.replics -c $gencode -l $params.scale_tree --write_anc
// echo "Mutation samples generated"

// for fasta_file in seqfile_sample-*.fasta
// do
// 	echo "Processing \$fasta_file"
// 	alignment2iqtree_states.py \$fasta_file  \${fasta_file}.state
// 	collect_mutations.py --tree ${tree}.ingroup --states  \${fasta_file}.state --gencode $gencode --syn $params.syn4f --no-mutspec --outdir mout --force
// 	cat mout/run.log >> pyvolve_${label}.log
// 	echo -e "\n\n">> pyvolve_${label}.log
// 	cat mout/mutations.tsv >  \${fasta_file}.mutations
// done
// echo "Mutations extraction done"

// concat_mutations.py seqfile_sample-*.fasta.mutations mutations_${label}_pyvolve.tsv
// echo "Mutations concatenation done"

// calculate_mutspec.py -b mutations_${label}_pyvolve.tsv -e mout/expected_freqs.tsv -o . \
// 	-l ${label}_simulated --exclude OUTGRP,ROOT --syn $params.syn4f $params.all --mnum192 $params.mnum192 --plot -x pdf \
// 	--substract12 \$spectra12 \$substractor192
// echo "Mutational spectrum calculated"

// """
// }


// process count_prior_mutations {

// publishDir params.outdir, overwrite: true, mode: 'copy',
// 	saveAs: {filename ->
// 	if (filename =~ /mutnumbers.tsv$/) "mutspec_tables/$filename"
// }

// input:
//  file seqs from g_433_multipleFasta_g_420

// output:
//  file "mutnumbers.tsv"  into g_420_outputFileTSV

// """
// /opt/dolphin/scripts/mutnumbers.pl $seqs 1>mutnumbers.tsv

// """
// }



workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
