if (!params.sequence){
	println "ERROR: Specify input nucleotide multifasta file"
	exit 1
}
if (!params.gencode){
	println "ERROR: Specify gencode number (e.g. 1,2,3 etc.)"
	exit 1
}
if (!params.species_name){
	println "ERROR: Specify species name"
	exit 1
}
if (!params.use_macse){params.use_macse = "false"} 
if (!params.verbose){params.verbose = "false"} 
if (!params.internal){params.internal = "false"} 
if (!params.terminal){params.terminal = "false"} 
if (!params.branch_spectra){params.branch_spectra = "false"}
if (!params.save_exp_mutations){params.save_exp_mutations = "false"}
if (!params.exclude_cons_sites){params.exclude_cons_sites = "false"}
if (!params.uncertainty_coef){params.uncertainty_coef = "false"}
if (!params.njobs){params.njobs = "1"}
THREADS = params.njobs

// if DB is NT, check that genus taxid specified
if (params.DB.endsWith("nt")){
	if (!params.genus_taxid){
		println "ERROR: Specify genus taxid when run pipeline on NT database"
		exit 1
	}
}

// TODO add specific params logs
if (params.verbose == 'true') {
	println ""
	println "PARAMETERS:"
	println "all: ${params.all}"
	println "syn: true"
	println "syn4f: ${params.syn4f}"
	println "Minimal number of mutations to save 192-component spectrum (mnum192): ${params.mnum192}"
	println "Use probabilities: ${params.use_probabilities}"
	if (params.use_probabilities == 'true'){
		println "Mutation probability cutoff: ${params.proba_cutoff}"
	}
	println "Run simulation: ${params.run_simulation}"
	if (params.run_simulation == 'true'){
		println "Replics in simulation: ${params.replics}"
		println "Tree scaling coefficient: ${params.scale_tree}"
	}
	println "Threads: ${THREADS}"
	println "Run treeshrink: ${params.run_shrinking}"
	if (params.run_shrinking == 'true'){
		println "Shrinking Quantile: ${params.quantile}"
	}
	println "Exclude conservative sites: ${params.exclude_cons_sites}"
	println "internal branches spectra: ${params.internal}"
	println "terminal branches spectra: ${params.terminal}"
	println "branch-scpecific spectra: ${params.branch_spectra}"
	println ""
}

params.all_arg   = params.all == "true" ? "--all" : ""
params.syn4f_arg = params.syn4f == "true" ? "--syn4f" : ""
params.proba_arg = params.use_probabilities == "true" ? "--proba" : ""

Channel.value(params.DB).set{g_15_commondb_path_g_406}
Channel.value(params.genus_taxid).set{genus_taxid_value}
query_protein_sequence = file(params.sequence, type: 'any') 
Channel.value(params.gencode).into{g_220_gencode_g_406;g_396_gencode_g_410;g_396_gencode_g_411;g_396_gencode_g_422;g_396_gencode_g_423;g_396_gencode_g_433}
Channel.value("false").set{aligned_param}
Channel.value(params.species_name).set{g_1_species_name_g_415}
Channel.value(params.use_macse).set{use_macse_param}


process query_qc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /input_seq_char_counts.log$/) "logs/$filename"
}

input:
 file query from query_protein_sequence

output:
 file "query_single.fasta"  into g_398_multipleFasta_g_406
 file "input_seq_char_counts.log" optional true  into g_398_logFile

"""
echo "Init query QC" >&2
if [ `grep -c ">" $query` -ne 1 ]; then
	echo "ERROR: Query fasta must contain single amino acid sequence" >&2
	exit 1
else
	echo "INFO: Number of sequences: `grep -c '>' $query`" >&2
fi

grep -v  ">" $query | grep -o . | sort | uniq -c | sort -nr > input_seq_char_counts.log
if [ `head -n 4 input_seq_char_counts.log | grep -Ec "[ACGT]"` -lt 4 ] || [ `grep -Ec "[EFILPQU]" input_seq_char_counts.log` -ne 0 ]; then
	echo "INFO: It's probably amino asid sequence" >&2
else
	echo "ERROR: Query fasta must contain single amino acid sequence" >&2
	exit 1
fi

mv $query query_single.fasta

"""
}


NSEQS_LIMIT=33000

process tblastn_and_seqs_extraction {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /sampled_sequences.fasta$/) "$filename"
	else if (filename =~ /report.blast$/) "logs/$filename"
	else if (filename =~ /.blast_output_*_filtered.tsv$/) "logs/$filename"
	else if (filename =~ /headers_mapping.txt$/) "$filename"
	else if (filename =~ /encoded_headers.txt$/) "$filename"
	else if (filename =~ /.*.taxids$/) "logs/$filename"
}

input:
 file query from g_398_multipleFasta_g_406
 val species_name from g_1_species_name_g_415
 val genus_taxid from genus_taxid_value
 val gencode from g_220_gencode_g_406
 val DB from g_15_commondb_path_g_406

output:
 file "sampled_sequences.fasta" into g_415_multipleFasta_g_409
 file "report.blast" optional true
 file "blast_output_*_filtered.tsv" optional true
 file "headers_mapping.txt" optional true
 file "encoded_headers.txt" optional true
 file "*.taxids" optional true

script:
"""
nseqs=500
outfmt="6 saccver pident length qlen gapopen sstart send evalue bitscore sframe"

if [[ $DB == *nt ]]; then
	echo "INFO: Collecting taxids" >&2
	if [[ "$species_name" == Homo* ]]; then
		get_species_taxids.sh -t 9604 > genus.taxids
	else
		get_species_taxids.sh -t $genus_taxid  > genus.taxids
	fi
	sleep 2
	get_species_taxids.sh -n "$species_name" | head -n 5 > sp_tax.info
	if [[ `grep "rank : species" sp_tax.info` ]]; then 
		raw_sp_taxid=`grep Taxid sp_tax.info`
		species_taxid="\${raw_sp_taxid#*Taxid : }"
	fi
	sleep 1
	get_species_taxids.sh -t \$species_taxid > species.taxids

	grep -v -f species.taxids genus.taxids > relatives.taxids

	echo "INFO: Checking number of taxids for outgroup" >&2
	if [ `wc -l relatives.taxids | cut -f 1 -d ' '` -eq 0 ]; then
		echo "ERROR: there are no taxids that can be used as outgroup." >&2
		echo "Maybe this species is single in the genus, so pipeline can build incorrect phylogeneti tree" >&2
		echo "due to potential incorrect tree rooting. You can select sequences and outgroup manually and" >&2
		echo "run pipeline on your nucleotide sequences" >&2
		exit 1
	fi

	echo "INFO: Blasting species sequences in the nt" >&2
	tblastn -db $DB -db_gencode $gencode -max_target_seqs \$nseqs \
			-query $query -out blast_output_species.tsv -evalue 0.00001 \
			-num_threads $THREADS -taxidlist species.taxids \
			-outfmt "\$outfmt"

	echo "INFO: Filtering out bad hits: ident < 80, coverage < 0.8, gap_opened > 0" >&2
	awk '\$5 == 0 && \$2 > 80 && \$3 / \$4 > 0.8' blast_output_species.tsv > blast_output_species_filtered.tsv

	echo "INFO: Checking required number of hits" >&2
	nhits=`wc -l blast_output_species_filtered.tsv | cut -f 1 -d ' '`
	if [ \$nhits -lt 10 ]; then
		echo "ERROR: there are only \$nhits valuable hits in the database for given query," >&2
		echo "but needed at least 10" >&2
		exit 1
	fi

	echo "INFO: Preparing coords for nucleotide sequences extraction" >&2
	# entry|range|strand, e.g. ID09 x-y minus
	awk '{print \$1, \$6, \$7, (\$NF ~ /^-/) ? "minus" : "plus"}' blast_output_species_filtered.tsv > raw_coords.txt
	awk '\$2 > \$3 {print \$1, \$3 "-" \$2, \$4}' raw_coords.txt > coords.txt
	awk '\$3 > \$2 {print \$1, \$2 "-" \$3, \$4}' raw_coords.txt >> coords.txt

	echo "INFO: Getting nucleotide sequences" >&2
	blastdbcmd -db $DB -entry_batch coords.txt -outfmt %f -out nucleotide_sequences.fasta

	echo "INFO: Checking required number of extracted seqs" >&2
	nseqs=`grep -c '>' nucleotide_sequences.fasta`
	if [ \$nseqs -lt 10 ]; then
		echo "ERROR: cannot extract more than \$nseqs seqs from the database for given query, but needed at least 10" >&2
		exit 1
	fi

	echo -e "INFO: Blasting for outgroup search" >&2
	tblastn -db $DB -db_gencode $gencode -max_target_seqs 10 \
			-query $query -out blast_output_genus.tsv -evalue 0.00001 \
			-num_threads $THREADS -taxidlist relatives.taxids \
			-outfmt "\$outfmt"
	
	echo "INFO: Filtering out bad hits: ident > 80, coverage > 0.8, gap_opened == 0" >&2
	awk '\$5 == 0 && \$2 > 80 && \$3 / \$4 > 0.8' blast_output_genus.tsv | sort -rk 9 > blast_output_genus_filtered.tsv

	echo "INFO: Checking required number of hits for outgroup" >&2
	if [ `wc -l blast_output_genus_filtered.tsv | cut -f 1 -d ' '` -eq 0 ]; then
		echo "ERROR: there are no hits in the database that could be used as outgroup" >&2
		exit 1
	fi

	echo "INFO: Preparing genus coords for nucleotide sequences extraction" >&2
	# entry|range|strand, e.g. ID09 x-y minus
	head -n 1 blast_output_genus_filtered.tsv | awk '{print \$1, \$6, \$7, (\$NF ~ /^-/) ? "minus" : "plus"}' > raw_coords_genus.txt
	awk '\$2 > \$3 {print \$1, \$3 "-" \$2, \$4}' raw_coords_genus.txt >  coords_genus.txt
	awk '\$3 > \$2 {print \$1, \$2 "-" \$3, \$4}' raw_coords_genus.txt >> coords_genus.txt

	echo "INFO: Getting nucleotide sequences of potential outgroups" >&2
	blastdbcmd -db $DB -entry_batch coords_genus.txt -outfmt %f -out outgroup_sequence.fasta
	
	echo "INFO: Sequences headers encoding" >&2
	ohead=`head -n 1 outgroup_sequence.fasta`
	cat outgroup_sequence.fasta nucleotide_sequences.fasta > sample.fasta
	multifasta_coding.py -a sample.fasta -g "\${ohead:1}" -o sampled_sequences.fasta -m encoded_headers.txt

else
	report=report.blast
	while true
	do   
		echo "INFO: Blasting in midori2 database; nseqs=\$nseqs" >&2
		tblastn -db $DB -db_gencode $gencode -num_descriptions \$nseqs -num_alignments \$nseqs \
				-query $query -out \$report -num_threads $THREADS
		
		if [ `grep -c "No hits found" \$report` -eq 0 ]; then 
			echo "INFO: some hits found in the database for given query" >&2
		else
			echo "ERROR: there are no hits in the database for given query" >&2
			exit 1
		fi

		if [ `grep -c "$species_name" \$report` -ge \$((nseqs * 2 - 10)) ]; then
			nseqs=\$((nseqs * 4))
			if [ \$nseqs -gt $NSEQS_LIMIT ]; then
				echo "UNEXPECTED ERROR: database cannot contain more than $NSEQS_LIMIT sequences of one gene of some species" >&2
				exit 1
			fi
			echo "INFO: run blasting again due to abcence of different species; nseqs=\$nseqs" >&2
		else
			echo "SUCCESS: other species for outgroup selection is in hits" >&2
			break
		fi
	done

	mview -in blast -out fasta \$report 1>raw_sequences.fasta

	/opt/scripts_latest/header_sel_mod3.pl raw_sequences.fasta "$species_name" 1>useless_seqs.fasta 2>headers_mapping.txt

	/opt/scripts_latest/nuc_coding_mod.pl headers_mapping.txt $DB 1>sampled_sequences.fasta

fi
"""
}


process duplicates_filtration {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /seqs_unique.fasta$/) "$filename"
	else if (filename =~ /report_(yes|no).log$/) "logs/$filename"
}

input:
 file seqs from g_415_multipleFasta_g_409

output:
 file "seqs_unique.fasta"  into g_409_multipleFasta_g418_428
 file "report_{yes,no}.log"  into g_499_logFile

"""
/opt/scripts_latest/codon_alig_unique.pl $seqs 1>seqs_unique.fasta
if [ -f report_no.txt ]; then
	mv report_no.txt report_no.log
	echo "ERROR: Cannot find in the database required number of sequences" >&2
	echo "(10) for given gene of the species!" >&2
	cat report_no.log >&2
	exit 1
else
	mv report_yes.txt report_yes.log
fi
"""
}


outgrp = "OUTGRP"
min_input_nseqs = 10

process nucleotide_fasta_qc {

input:
 file query from g_409_multipleFasta_g418_428

output:
 file "sequences.fasta"  into g_428_multipleFasta_g_433
 file "char_numbers.log"

"""
if [ `grep -c ">" $query` -lt $min_input_nseqs ]; then
	echo "Number of sequences must be >= $min_input_nseqs" >&2
	exit 1
fi

if [ `grep -E -c ">$outgrp" $query` -ne 1 ]; then
	echo "Cannot find outgroup header in the alignment." >&2
	echo "Probably outgroup sequence for this gene cannot be found in the database." >&2
	exit 1
fi

grep -v  ">" $query | grep -o . | sort | uniq -c | sort -nr > char_numbers.log
if [ `head -n 5 char_numbers.log | grep -Ec "[ACGTacgt]"` -ge 3 ] && [ `grep -Ec "[EFILPQU]" char_numbers.log` -eq 0 ]; then
	echo "All right" >&2
else
	echo "Query fasta must contain nucleotides" >&2
	exit 1
fi

# this script do nothing, TODO drop it. This version of pipeline is for single protein!
multifasta_coding.py -a $query -g $outgrp -o sequences.fasta -m nothing.txt
"""
}


process MSA {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /msa_nuc.fasta$/) "$filename"
}

input:
 val aligned from aligned_param
 file seqs from g_428_multipleFasta_g_433
 val gencode from g_396_gencode_g_433
 val use_macse from use_macse_param

output:
 file "msa_nuc.fasta"  into g_433_multipleFasta_g_420, g_433_multipleFasta_g_421, g_433_multipleFasta_g_422, g_433_multipleFasta_g_424, g_433_multipleFasta_g_425
 file "aln_aa.fasta" optional true  into g_433_multipleFasta_lol

"""
if [ $aligned = true ]; then
	# TODO add verification of aligment and delete this useless argument!!!!
	echo "No need to align!" >&2
	cp $seqs aln.fasta
elif [ $aligned = false ]; then
	if [ $use_macse = true ]; then
		echo "Use macse as aligner" >&2
		java -jar /opt/macse_v2.06.jar -prog alignSequences -gc_def $gencode \
			-out_AA aln_aa.fasta -out_NT aln.fasta -seq $seqs
	else
		echo "Use mafft as aligner" >&2
		# NT2AA
		java -jar /opt/macse_v2.06.jar -prog translateNT2AA -seq $seqs \
			-gc_def $gencode -out_AA translated.faa
		#ALN AA
		mafft --thread $THREADS translated.faa > translated_aln.faa
		#AA_ALN --> NT_ALN
		java -jar /opt/macse_v2.06.jar -prog reportGapsAA2NT \
			-align_AA translated_aln.faa -seq $seqs -out_NT aln.fasta

	fi
else
	echo "ERROR: Inappropriate value for 'aligned' parameter; must be 'true' or 'false'" >&2
	exit 1
fi

echo "Do quality control" >&2
/opt/scripts_latest/macse2.pl aln.fasta msa_nuc.fasta
"""
}


process write_readme {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /readme.txt$/) "$filename"
}

input:
 file something from g_433_multipleFasta_g_425
output:
 file "readme.txt"

"""
cat > readme.txt <<- EOM
Output structure:

.
├── final_tree.nwk						# Final phylogenetic tree
├── seqs_unique.fasta					# Filtered orthologous sequences
├── msa_nuc.fasta						# Verified multiple sequence alignment
├── headers_mapping.txt					# Encoded headers of sequences
├── encoded_headers.txt					# Encoded headers of sequences (v2 for different versions of input)
├── logs/
│   ├── char_numbers.log				# Character composition of input sequence (used in query QC)
│   ├── report.blast					# Tblastn output during orthologs search
│   ├── *.taxids						# Taxids used in taxa-specific blasing in nt; relatives.taxids contains 
│	│									# 	other species from the genus of query and used for outgroup selection
│   ├── iqtree.log						# IQ-TREE logs during phylogenetic tree inference
│   ├── iqtree_report.log				# IQ-TREE report during phylogenetic tree inference
│   ├── iqtree_treeshrink.log			# TreeShrink logs
│   ├── iqtree_pruned_nodes.log			# Nodes pruned from tree by TreeShrink
│   ├── iqtree_anc.log					# IQ-TREE logs during ancestral reconstrution
│   ├── iqtree_anc_report.log			# IQ-TREE report during ancestral reconstrution
│   ├── iqtree_mut_extraction.log		# Logs during mutation extraction process
│   └── branches.txt					# Tree branch lenghts
├── figures
│   ├── ms12syn.pdf						# Barplot with  12-component spectrum on synonymous mutations
│   └── ms192syn.pdf					# Barplot with 192-component spectrum on synonymous mutations
├── tables
│   ├── rates.tsv						# Site rates categories for an alignment
│   ├── expected_freqs.tsv				# Frequencies of substitutions for each tree node genome
│   ├── mean_expexted_mutations.tsv		# Averaged frequencies of substitutions for entire tree
│   ├── ms12syn.tsv						# table with 12-component spectrum on synonymous mutations
│   ├── ms192syn.tsv					# table with 192-component spectrum on synonymous mutations
│   └── observed_mutations.tsv			# Recontructed mutations
EOM
"""
}


process fasta2states_table {

input:
 file aln from g_433_multipleFasta_g_421

output:
 file "leaves_states.state"  into g_421_state_g_410

"""
alignment2iqtree_states.py $aln leaves_states.state
"""
}


process ML_tree_inference {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 file mulal from g_433_multipleFasta_g_420

output:
 set val("iqtree"), file("iqtree.nwk")  into g_409_tree_g_315
 file "*.log" optional true  into g_409_logFile

errorStrategy 'retry'
maxRetries 3

script:
"""
iqtree2 -s $mulal -m $params.iqtree_model -nt $THREADS --prefix phylo
mv phylo.treefile iqtree.nwk
mv phylo.iqtree iqtree_report.log
mv phylo.log iqtree.log
"""
}


process shrink_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 set val(name), file(tree) from g_409_tree_g_315

output:
 set val("${name}_shrinked"), file("${name}_shrinked.nwk")  into g_315_tree_g_302, g_315_tree_g_132
 file "*.log" optional true into g_315_logFile

"""
if [ $params.run_shrinking = true ] && [ `nw_stats $tree | grep leaves | cut -f 2` -gt 8 ]; then
	run_treeshrink.py -t $tree -O treeshrink -o . -q $params.quantile -x OUTGRP
	mv treeshrink.nwk ${name}_shrinked.nwk
	mv treeshrink_summary.txt ${name}_treeshrink.log
	mv treeshrink.txt ${name}_pruned_nodes.log
else
	cat $tree > ${name}_shrinked.nwk
	if [ $params.run_shrinking = true ]; then
		echo "Shrinking are useless on such a small number of sequences" > ${name}_treeshrink.log
	fi
fi
"""
}


process terminal_branch_lengths_qc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /branches.txt$/) "logs/$filename"
}

input:
 set val(name), file(tree) from g_315_tree_g_132

output:
 set val(name), file("branches.txt")  into g_132_branches

"""
nw_distance -m p -s f -n $tree | sort -grk 2 1> branches.txt

if [ `grep OUTGRP branches.txt | cut -f 2 | python3 -c "import sys; print(float(sys.stdin.readline().strip()) > 0)"` = False ]; then
	cat "branches.txt"
	echo "Something went wrong: outgroup is not furthest leaf in the tree" >&2
	exit 1
fi
"""
}


process tree_rooting_iqtree {

input:
 set val(name), file(tree) from g_315_tree_g_302

output:
 set val("${name}_rooted"), file("*.nwk")  into g_302_tree_g_326

"""
nw_reroot -l $tree OUTGRP 1>${name}_rooted.nwk
"""
}


estimate_rates = params.exclude_cons_sites == "true" ? "--rate" : ""

process ASR {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /final_tree.nwk$/) "$filename"
	else if (filename =~ /rates.tsv$/) "tables/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
}

input:
 file mulal from g_433_multipleFasta_g_424
 set val(namet), file(tree) from g_302_tree_g_326

output:
 set val("iqtree"), file("final_tree.nwk")  into g_326_tree_g_410, g_326_tree_g_422
 set val("iqtree"), file("iqtree_anc.state")  into g_326_state_g_410
 path "rates.tsv" into g_326_ratefile
 file "*.log"  into g_326_logFile

errorStrategy 'retry'
maxRetries 3

script:
"""
nw_labels -I $tree | sed 's/\$/\$/' > leaves.txt
select_records.py -f leaves.txt --fmt fasta $mulal msa_filtered.fasta

iqtree2 -te $tree -s msa_filtered.fasta -m $params.iqtree_anc_model -asr -nt $THREADS --prefix anc $estimate_rates
if [ ! -f anc.rate ]; then
	touch anc.rate
fi
mv anc.rate rates.tsv
mv anc.iqtree iqtree_anc_report.log
mv anc.log iqtree_anc.log
nw_reroot anc.treefile OUTGRP | sed 's/;/ROOT;/' > final_tree.nwk

iqtree_states_add_part.py anc.state iqtree_anc.state
"""
}


save_exp_muts = params.save_exp_mutations == "true" ? "--save-exp-muts" : ""
use_uncertainty_coef = params.uncertainty_coef == "true" ? "--phylocoef" : "--no-phylocoef"

process mutations_extraction {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "tables/$filename"
	else if (filename =~ /.*.log$/) "logs/$filename"
	else if (filename =~ /.*.pdf$/) "figures/$filename"
}

input:
 set val(namet), file(tree) from g_326_tree_g_410
 set val(label), file(states1) from g_326_state_g_410
 path rates from g_326_ratefile
 file states2 from g_421_state_g_410
 val gencode from g_396_gencode_g_410

output:
 file "*.tsv"  into g_410_outputFileTSV
 file "*.log"  into g_410_logFile
 file "*.pdf"  into g_410_outputFilePdf
 file "ms*syn.txt" optional true  into g_410_outputFileTxt_g_422

"""
if [ $params.exclude_cons_sites = true ]; then 
	collect_mutations.py --tree $tree --states $states1 --states $states2 \
		--gencode $gencode --syn $params.syn4f_arg $params.proba_arg --no-mutspec \
		--pcutoff $params.proba_cutoff --mnum192 $params.mnum192 \
		--rates $rates --cat-cutoff $params.cons_cat_cutoff \
		--outdir mout $save_exp_muts $use_uncertainty_coef
else
	collect_mutations.py --tree $tree --states $states1 --states $states2 \
		--gencode $gencode --syn $params.syn4f_arg $params.proba_arg --no-mutspec \
		--pcutoff $params.proba_cutoff --mnum192 $params.mnum192 \
		--outdir mout $save_exp_muts $use_uncertainty_coef
fi
mv mout/* .
mv mutations.tsv observed_mutations.tsv
mv run.log mut_extraction.log

calculate_mutspec.py -b observed_mutations.tsv -e expected_freqs.tsv -o . \
	--exclude OUTGRP,ROOT --mnum192 $params.mnum192 $params.proba_arg \
	--proba_min $params.proba_cutoff --syn $params.syn4f_arg $params.all_arg --plot -x pdf

if [ $params.internal = true ]; then
	calculate_mutspec.py -b observed_mutations.tsv -e expected_freqs.tsv -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 $params.proba_arg \
		--proba_min $params.proba_cutoff --syn $params.syn4f_arg $params.all_arg --plot -x pdf --subset internal
	rm mean_expexted_mutations_internal.tsv
fi
if [ $params.terminal = true ]; then
	calculate_mutspec.py -b observed_mutations.tsv -e expected_freqs.tsv -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 $params.proba_arg \
		--proba_min $params.proba_cutoff --syn $params.syn4f_arg $params.all_arg --plot -x pdf --subset terminal
	rm mean_expexted_mutations_terminal.tsv
fi
if [ $params.branch_spectra = true ]; then
	calculate_mutspec.py -b observed_mutations.tsv -e expected_freqs.tsv -o . \
        --exclude OUTGRP,ROOT --mnum192 $params.mnum192 $params.proba_arg \
		--proba_min $params.proba_cutoff --syn $params.syn4f_arg $params.all_arg --branches
fi

cp ms12syn.tsv ms12syn.txt
if [-f ms192syn.tsv ]; then
	cp ms192syn.tsv ms192syn.txt
fi
"""
}


if (params.verbose == 'true') {
	workflow.onComplete {
	println "##Completed at: ${workflow.complete}; Duration: ${workflow.duration}; Success: ${workflow.success ? 'OK' : 'failed' }"
	}
}
