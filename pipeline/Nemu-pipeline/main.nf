THREADS = 4

if (!params.species_name){params.species_name = ""} 
if (!params.sequence){params.sequence = ""} 
if (!params.Mt_DB){params.Mt_DB = ""} 
if (!params.gencode){params.gencode = ""} 

Channel.value(params.species_name).set{g_1_species_name_g_414}
g_2_multipleFasta_g_398 = file(params.sequence, type: 'any') 
Channel.value(params.Mt_DB).into{g_15_commondb_path_g_406;g_15_commondb_path_g_415}
Channel.value(params.gencode).into{g_220_gencode_g_406;g_220_gencode_g418_410;g_220_gencode_g418_411;g_220_gencode_g418_422;g_220_gencode_g418_423;g_220_gencode_g418_433}


process query_qc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /char_numbers.log$/) "logs/$filename"
}

input:
 file query from g_2_multipleFasta_g_398

output:
 file "query_single.fasta"  into g_398_multipleFasta_g_406
 file "char_numbers.log" optional true  into g_398_logFile

"""
if [ `grep -c ">" $query` -ne 1 ]; then
	echo "Query fasta must contain single amino acid sequence"
	exit 1
fi

grep -v  ">" $query | grep -o . | sort | uniq -c | sort -nr > char_numbers.log
if [ `head -n 4 char_numbers.log | grep -Ec "[ACGT]"` -lt 4 ] || [ `grep -Ec "[EFILPQU]" char_numbers.log` -ne 0 ]; then
	echo "It's probably amino asid sequence"
else
	echo "Query fasta must contain single amino acid sequence"
	exit 1
fi

mv $query query_single.fasta

"""
}

nseqs = params.tblastn.nseqs


process tblastn {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /report.blast$/) "logs/$filename"
	else if (filename =~ /no_hits.log$/) "logs/$filename"
}

input:
 val gencode from g_220_gencode_g_406
 file query from g_398_multipleFasta_g_406
 val DB from g_15_commondb_path_g_406

output:
 file "report.blast"  into g_406_blast_output_g_412
 file "no_hits.log" optional true  into g_406_logFile

script:
"""
tblastn -db $DB -db_gencode $gencode -num_descriptions $nseqs -num_alignments $nseqs -query $query -out report.blast -num_threads $THREADS

if [ `grep -c "No hits found" report.blast` -eq 0 ]; then 
	echo "Found hits in the database for given query"
else
	echo "There are no hits in the database for given query" > no_hits.log
	cat no_hits.log
	exit 1
fi
"""

}


process mview {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /raw_sequences.fasta$/) "tmp/$filename"
}

input:
 file blast_report from g_406_blast_output_g_412

output:
 file "raw_sequences.fasta"  into g_412_multipleFasta_g_414

"""
mview -in blast -out fasta $blast_report 1>raw_sequences.fasta

"""
}


process extract_outgroup {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /useless_seqs.fasta$/) "tmp/$filename"
	else if (filename =~ /headers_mapping.txt$/) "sequences/$filename"
}

input:
 val SPNAME from g_1_species_name_g_414
 file query_out_fasta from g_412_multipleFasta_g_414

output:
 file "useless_seqs.fasta"  into g_414_multipleFasta
 file "headers_mapping.txt"  into g_414_outputFileTxt_g_415

"""
/opt/dolphin/scripts/header_sel_mod3.pl $query_out_fasta "$SPNAME" 1>useless_seqs.fasta 2>headers_mapping.txt

"""
}


process extract_sequences {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /sequences.fasta$/) "tmp/$filename"
}

input:
 file hash from g_414_outputFileTxt_g_415
 val DB from g_15_commondb_path_g_415

output:
 file "sequences.fasta"  into g_415_multipleFasta_g_409

"""
/opt/dolphin/scripts/nuc_coding_mod.pl $hash $DB 1>sequences.fasta

"""
}


process drop_dublicates {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /seqs_unique.fasta$/) "sequences/$filename"
	else if (filename =~ /report_(yes|no).log$/) "logs/$filename"
}

input:
 file seqs from g_415_multipleFasta_g_409

output:
 file "seqs_unique.fasta"  into g_409_multipleFasta_g418_428
 file "report_{yes,no}.log"  into g_409_logFile

"""
/opt/dolphin/scripts/codon_alig_unique.pl $seqs 1>seqs_unique.fasta
if [ -f report_no.txt ]; then
	mv report_no.txt report_no.log
	echo "ERROR! Cannot find in the database required number of sequsences (10) for given species!"
	cat report_no.log
	exit 1
else
	mv report_yes.txt report_yes.log
fi

"""
}


process pass_variables {


output:
 val "OUTGRP"  into g_420_outgroup_g418_428
 val "single"  into g_420_mode_g418_410, g_420_mode_g418_411, g_420_mode_g418_422, g_420_mode_g418_423
 val "false"  into g_420_type_g418_433

"""
echo Nothing
"""
}


process Nemu_tail_fasta_qc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /species_mapping.txt$/) "sequences/$filename"
}

input:
 file query from g_409_multipleFasta_g418_428
 val outgrp from g_420_outgroup_g418_428

output:
 file "sequences.fasta"  into g418_428_multipleFasta_g418_433
 file "species_mapping.txt" optional true  into g418_428_outputFileTxt
 file "char_numbers.log"  into g418_428_logFile

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


process Nemu_tail_macse {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /alignment_checked.fasta$/) "sequences/$filename"
}

input:
 val aligned from g_420_type_g418_433
 file seqs from g418_428_multipleFasta_g418_433
 val gencode from g_220_gencode_g418_433

output:
 file "alignment_checked.fasta"  into g418_433_multipleFasta_g418_420, g418_433_multipleFasta_g418_421, g418_433_multipleFasta_g418_423, g418_433_multipleFasta_g418_422, g418_433_multipleFasta_g418_424
 file "aln_aa.fasta" optional true  into g418_433_multipleFasta

"""
if [ $aligned == "true" ]; then
	echo "No need to align!"
	cp $seqs aln.fasta
elif [ $aligned == "false" ]; then
	if [ `grep -c ">" $seqs` -le 1000 ]; then
		java -jar /opt/macse_v2.06.jar -prog alignSequences -gc_def $gencode -out_AA aln_aa.fasta -out_NT aln.fasta -seq $seqs
	else
		echo "Too many seuqences (>1000), need to use mafft as aligner"
		# NT2AA
		java -jar /opt/macse_v2.06.jar -prog translateNT2AA -seq $seqs -gc_def $gencode -out_AA translated.faa
		#ALN AA
		mafft --thread $THREADS translated.faa > translated_aln.faa
		#AA_ALN --> NT_ALN
		java -jar /opt/macse_v2.06.jar -prog reportGapsAA2NT -align_AA translated_aln.faa -seq $seqs -out_NT aln.fasta
	fi
else
	echo "Inappropriate values for 'aligned'"
	exit 1
fi

echo "Do quality control"
/opt/dolphin/scripts/macse2.pl aln.fasta alignment_checked.fasta

"""
}


process Nemu_tail_convert_alignment_to_phylip {

input:
 file aln from g418_433_multipleFasta_g418_424

output:
 set val("aln"), file("aln.phy")  into g418_424_phylip_g418_130, g418_424_phylip_g418_326, g418_424_phylip_g418_409, g418_424_phylip_g418_430

"""
java -jar /opt/readseq.jar -a -f Phylip -o aln.phy $aln

"""
}

run_IQTREE = params.Nemu_tail_iqtree_build_tree.run_IQTREE
iqtree_model = params.Nemu_tail_iqtree_build_tree.iqtree_model
quantile = params.Nemu_tail_iqtree_build_tree.quantile
//* @style @condition:{run_IQTREE="true", iqtree_model, quantile}

params.IQTREE_model = iqtree_model
params.quantile = quantile

//--mset 6.7a,6.7b,6.8a,8.8,8.10a,8.16,8.17,8.18,10.12,10.34,12.12 -mrate FO+R6+I

process Nemu_tail_iqtree_build_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 set val(name), file(mulal) from g418_424_phylip_g418_409

output:
 set val("iqtree"), file("iqtree.nwk")  into g418_409_tree_g418_315
 file "*.log"  into g418_409_logFile

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


process Nemu_tail_shrink_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 set val(name), file(tree) from g418_409_tree_g418_315

output:
 set val("${name}_shrinked"), file("${name}_shrinked.nwk")  into g418_315_tree_g418_302, g418_315_tree_g418_132
 file "*.log"  into g418_315_logFile

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


process Nemu_tail_extract_terminal_branch_lengths_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.branches$/) "IQTREE/$filename"
}

input:
 set val(name), file(tree) from g418_315_tree_g418_132

output:
 set val(name), file("${name}.branches")  into g418_132_branches

"""
nw_distance -m p -s l -n $tree | sort -grk 2 1> ${name}.branches

if [ `grep OUTGRP ${name}.branches | cut -f 2 | python3 -c "import sys; print(float(sys.stdin.readline().strip()) > 0)"` == False ]; then
	cat "${name}.branches"
	echo "Something went wrong: outgroup is not furthest leaf in the tree"
	exit 1
fi

"""
}


process Nemu_tail_rooting_iqtree_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.nwk$/) "tmp/$filename"
}

input:
 set val(name), file(tree) from g418_315_tree_g418_302

output:
 set val("${name}_rooted"), file("*.nwk")  into g418_302_tree_g418_326

"""
nw_reroot -l $tree OUTGRP 1>${name}_rooted.nwk

"""
}


process Nemu_tail_iqtree_anc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /iqtree_anc_tree.nwk$/) "IQTREE/$filename"
	else if (filename =~ /iqtree_anc.state$/) "IQTREE/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 set val(name), file(mulal) from g418_424_phylip_g418_326
 set val(namet), file(tree) from g418_302_tree_g418_326

output:
 set val("iqtree"), file("iqtree_anc_tree.nwk")  into g418_326_tree_g418_410, g418_326_tree_g418_422
 set val("iqtree"), file("iqtree_anc.state")  into g418_326_state_g418_410
 file "*.log"  into g418_326_logFile

errorStrategy 'retry'
maxRetries 3

script:
"""
nw_labels -I $tree > leaves.txt
filter_aln.py -a $mulal -l leaves.txt -o filtered_aln.phy

iqtree2 -te $tree -s filtered_aln.phy -m $params.IQTREE_model -asr -nt $THREADS --prefix anc
mv anc.iqtree iqtree_anc_report.log
# mv anc.state iqtree_anc.state
mv anc.log iqtree_anc.log
nw_reroot anc.treefile OUTGRP | sed 's/;/ROOT;/' > iqtree_anc_tree.nwk

iqtree_states_add_part.py anc.state iqtree_anc.state
"""

}

run_RAXML = params.Nemu_tail_raxml_build_tree.run_RAXML
raxml_model = params.Nemu_tail_raxml_build_tree.raxml_model
//* @style @condition:{run_RAXML="true", raxml_model}

params.RAxML_model = raxml_model


process Nemu_tail_raxml_build_tree {

input:
 set val(name), file(mulal) from g418_424_phylip_g418_130

output:
 set val("raxml"), file("raxml.nwk")  into g418_130_tree_g418_317

when:
run_RAXML == "true"

script:
"""
/opt/standard-RAxML-8.2.12/raxmlHPC-PTHREADS-SSE3 -x 987654 -p 987654 -T $THREADS -N 50 -f a -s $mulal -n Rax_tree -m $raxml_model
mv RAxML_bestTree.Rax_tree raxml.nwk
"""

}


process Nemu_tail_shrink_raxml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 set val(name), file(tree) from g418_130_tree_g418_317

output:
 set val("${name}_shrinked"), file("${name}_shrinked.nwk")  into g418_317_tree_g418_301, g418_317_tree_g418_109
 file "*.log"  into g418_317_logFile

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


process Nemu_tail_extract_terminal_branch_lengths_raxml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.branches$/) "RAxML/$filename"
}

input:
 set val(name), file(tree) from g418_317_tree_g418_109

output:
 set val(name), file("${name}.branches")  into g418_109_branches

"""
nw_distance -m p -s l -n $tree | sort -grk 2 1> ${name}.branches

if [ `grep OUTGRP ${name}.branches | cut -f 2 | python3 -c "import sys; print(float(sys.stdin.readline().strip()) > 0)"` == False ]; then
	cat "${name}.branches"
	echo "Something went wrong: outgroup is not furthest leaf in the tree"
	exit 1
fi

"""
}


process Nemu_tail_rooting_raxml_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.nwk$/) "tmp/$filename"
}

input:
 set val(name), file(tree) from g418_317_tree_g418_301

output:
 set val("${name}_rooted"), file("*.nwk")  into g418_301_tree_g418_430

"""
nw_reroot -l $tree OUTGRP 1>${name}_rooted.nwk

"""
}


process Nemu_tail_raxml_anc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /RAxML_nodeLabelledRootedTree.nwk$/) "RAxML/$filename"
	else if (filename =~ /RAxML_marginalAncestralProbabilities.state$/) "RAxML/$filename"
	else if (filename =~ /RAxML_anc_rec.log$/) "logs/$filename"
	else if (filename =~ /RAxML_marginalAncestralProbabilities.ANCESTORS.txt$/) "tmp/$filename"
}

input:
 set val(namet), file(tree) from g418_301_tree_g418_430
 set val(name), file(mulal) from g418_424_phylip_g418_430

output:
 set val("raxml"),file("RAxML_nodeLabelledRootedTree.nwk")  into g418_430_tree_g418_411, g418_430_tree_g418_423
 set val("raxml"), file("RAxML_marginalAncestralProbabilities.state")  into g418_430_state_g418_411
 file "RAxML_marginalAncestralStates.fasta"  into g418_430_multipleFasta
 file "RAxML_anc_rec.log"  into g418_430_logFile
 file "RAxML_marginalAncestralProbabilities.ANCESTORS.txt" optional true  into g418_430_outputFileTxt

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


mv RAxML_marginalAncestralProbabilities.ANCESTORS RAxML_marginalAncestralProbabilities.ANCESTORS.txt
"""
}


process Nemu_tail_terminal_genomes_states {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /leaves_states.state$/) "tmp/$filename"
}

input:
 file aln from g418_433_multipleFasta_g418_421

output:
 set val("leaves_states"), file("leaves_states.state")  into g418_421_state_g418_410, g418_421_state_g418_411

"""
alignment2iqtree_states.py $aln leaves_states.state
"""
}


process Nemu_tail_mutations_raxml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "tables/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
	else if (filename =~ /.*.pdf$/) "images/$filename"
}

input:
 set val(namet), file(tree) from g418_430_tree_g418_411
 set val(label), file(states1) from g418_430_state_g418_411
 set val(names2), file(states2) from g418_421_state_g418_411
 val nspecies from g_420_mode_g418_411
 val gencode from g_220_gencode_g418_411

output:
 file "*.tsv"  into g418_411_outputFileTSV
 file "*.log"  into g418_411_logFile
 file "*.pdf"  into g418_411_outputFilePdf
 file "ms*syn_${label}.txt" optional true  into g418_411_outputFileTxt_g418_423

"""
if [ $nspecies == "single" ]; then
    collect_mutations.py --tree $tree --states $states1 --states $states2 \
        --gencode $gencode --syn $params.syn4f $params.proba --no-mutspec --outdir mout
	
	mv mout/* .
    mv mutations.tsv observed_mutations_${label}.tsv
    mv run.log ${label}_mut_extraction.log

    calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_freqs.tsv -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 -l $label $params.proba --proba_min $params.proba_min --syn $params.syn4f $params.all --plot -x pdf
    calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_freqs.tsv -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 -l $label $params.proba --proba_min $params.proba_min --syn $params.syn4f $params.all --plot -x pdf --subset internal

    cp ms12syn_${label}.tsv ms12syn_${label}.txt
    if [-f ms192syn_${label}.tsv ]; then
	    cp ms192syn_${label}.tsv ms192syn_${label}.txt
	fi
	
elif [ $nspecies == "multiple" ]; then
    collect_mutations.py --tree $tree --states $states1 --states $states2 \
        --gencode $gencode --syn $params.syn4f $params.proba --outdir mout
        
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


process Nemu_tail_pyvolve_raxml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "tables/$filename"
	else if (filename =~ /.*.pdf$/) "images/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 file spectra from g418_411_outputFileTxt_g418_423
 set val(label), file(tree) from g418_430_tree_g418_423
 file mulal from g418_433_multipleFasta_g418_423
 val gencode from g_220_gencode_g418_423
 val nspecies from g_420_mode_g418_423

output:
 file "*.tsv"  into g418_423_outputFileTSV
 file "*.pdf"  into g418_423_outputFilePdf
 file "*.log"  into g418_423_logFile

when:
params.run_pyvolve == "true" && nspecies == "single"

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


process Nemu_tail_mutations_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "tables/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
	else if (filename =~ /.*.pdf$/) "images/$filename"
}

input:
 set val(namet), file(tree) from g418_326_tree_g418_410
 set val(label), file(states1) from g418_326_state_g418_410
 set val(names2), file(states2) from g418_421_state_g418_410
 val nspecies from g_420_mode_g418_410
 val gencode from g_220_gencode_g418_410

output:
 file "*.tsv"  into g418_410_outputFileTSV
 file "*.log"  into g418_410_logFile
 file "*.pdf"  into g418_410_outputFilePdf
 file "ms*syn_${label}.txt" optional true  into g418_410_outputFileTxt_g418_422

"""
if [ $nspecies == "single" ]; then
    collect_mutations.py --tree $tree --states $states1 --states $states2 \
        --gencode $gencode --syn $params.syn4f $params.proba --no-mutspec --outdir mout
	
	mv mout/* .
    mv mutations.tsv observed_mutations_${label}.tsv
    mv run.log ${label}_mut_extraction.log

    calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_freqs.tsv -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 -l $label $params.proba --proba_min $params.proba_min --syn $params.syn4f $params.all --plot -x pdf
    calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_freqs.tsv -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 -l $label $params.proba --proba_min $params.proba_min --syn $params.syn4f $params.all --plot -x pdf --subset internal

    cp ms12syn_${label}.tsv ms12syn_${label}.txt
    if [-f ms192syn_${label}.tsv ]; then
	    cp ms192syn_${label}.tsv ms192syn_${label}.txt
	fi
	
elif [ $nspecies == "multiple" ]; then
    collect_mutations.py --tree $tree --states $states1 --states $states2 \
        --gencode $gencode --syn $params.syn4f $params.proba --outdir mout
        
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


process Nemu_tail_pyvolve_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "tables/$filename"
	else if (filename =~ /.*.pdf$/) "images/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 file spectra from g418_410_outputFileTxt_g418_422
 set val(label), file(tree) from g418_326_tree_g418_422
 file mulal from g418_433_multipleFasta_g418_422
 val gencode from g_220_gencode_g418_422
 val nspecies from g_420_mode_g418_422

output:
 file "*.tsv"  into g418_422_outputFileTSV
 file "*.pdf"  into g418_422_outputFilePdf
 file "*.log"  into g418_422_logFile

when:
params.run_pyvolve == "true" && nspecies == "single"

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


process Nemu_tail_count_prior_mutations {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /mutnumbers.tsv$/) "tables/$filename"
}

input:
 file seqs from g418_433_multipleFasta_g418_420

output:
 file "mutnumbers.tsv"  into g418_420_outputFileTSV

"""
/opt/dolphin/scripts/mutnumbers.pl $seqs 1>mutnumbers.tsv

"""
}


process Nemu_tail_mut_processing_params {



syn4f = params.Nemu_tail_mut_processing_params.syn4f
all = params.Nemu_tail_mut_processing_params.all
mnum192 = params.Nemu_tail_mut_processing_params.mnum192
use_probabilities = params.Nemu_tail_mut_processing_params.use_probabilities
proba_min = params.Nemu_tail_mut_processing_params.proba_min
simulation = params.Nemu_tail_mut_processing_params.simulation
replics = params.Nemu_tail_mut_processing_params.replics
scale_tree = params.Nemu_tail_mut_processing_params.scale_tree
site_rates = params.Nemu_tail_mut_processing_params.site_rates
categories = params.Nemu_tail_mut_processing_params.categories
cat_cutoff = params.Nemu_tail_mut_processing_params.cat_cutoff

//* @style @multicolumn:{syn4f, all, mnum192}, {use_probabilities, proba_min}, {simulation, replics, scale_tree}, {site_rates, categories, cat_cutoff} @condition:{simulation="true", replics, scale_tree}, {site_rates="true", categories, cat_cutoff}

println ""
println "Arguments:"
println "syn4f: ${syn4f}"
println "all: ${all}"
println "mnum192: ${mnum192}"
println "use_probabilities: ${use_probabilities}"
println "proba_min: ${proba_min}"
println "simulation: ${simulation}"
println "replics: ${replics}"
println "scale_tree: ${scale_tree}"
println "site_rates: ${site_rates}"
println "categories: ${categories}"
println "cat_cutoff: ${cat_cutoff}"
println "threads: ${THREADS}"
println ""

params.syn4f = syn4f == "true" ? "--syn4f" : ""
params.all = all == "true" ? "--all" : ""
params.mnum192 = mnum192
params.proba = use_probabilities == "true" ? "--proba" : ""
params.proba_min = proba_min
params.run_pyvolve = simulation
params.replics = replics
params.scale_tree = scale_tree
params.site_rates = site_rates
params.categories = categories
params.cat_cutoff = cat_cutoff

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
