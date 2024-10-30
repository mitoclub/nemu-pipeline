# syn4f="false" #//* @checkbox @description:"Run extraction of synonymous fourfold mutational spectra"
# use_probabilities="true" #//* @checkbox @description:"Use probabilities of nucleotides in mutational spectra calculation"
mnum192="16" #//* @input @description:"Number of mutation types (max 192) required to calculate and plot 192-component mutational spectra"
proba_cutoff="0.1" #//* @input @description:"Mutation probability cutoff: mutations with lower probability will not be considered in spectra calculation"
#//* @style @multicolumn:{syn4f, use_probabilities}, {mnum192, proba_cutoff}

# input
# val(namet), file(tree)
# val(label), file(states1)
# val(names2), file(states2)
label=iqtree
gencode=2
nspecies=single # single multiple

outdir=data/exposure/mammals_nd1/test_$nspecies
mkdir -p $outdir
cd $outdir

tree=../iqtree_anc_tree.nwk
states1=../iqtree_anc.state
states2=../leaves_states.state

# output
# "*.tsv"
# "*.log"
# "*.pdf"
# "ms12syn_${label}.txt"


if [ $nspecies == "single" ]; then
    # 3.collect_mutations.py --tree $tree --states $states1 --states $states2 \
    #     --gencode $gencode --syn --proba --no-mutspec --outdir .

    # mv mutations.tsv observed_mutations_${label}.tsv
    # mv expected_mutations.tsv expected_mutations.txt
    # mv run.log ${label}_mut_extraction.log

    calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_mutations.txt -o . \
        --exclude OUTGRP,ROOT --mnum192 $mnum192 -l $label --proba --proba_cutoff $proba_cutoff --syn --plot -x pdf
    calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_mutations.txt -o . \
        --exclude OUTGRP,ROOT --mnum192 $mnum192 -l $label --proba --proba_cutoff $proba_cutoff --syn --plot -x pdf --subset internal

    cp ms12syn_${label}.tsv ms12syn_${label}.txt

elif [ $nspecies == "multiple" ]; then
    3.collect_mutations.py --tree $tree --states $states1 --states $states2 \
        --gencode $gencode --syn --proba --outdir .

    mv mutations.tsv observed_mutations_${label}.tsv
    mv expected_mutations.tsv expected_mutations.txt
    mv run.log ${label}_mut_extraction.log

    python ~/mutspec-utils/scripts/plot_spectra.py -s mutspec12.tsv  -o ${label}_ms12.pdf  -t
    python ~/mutspec-utils/scripts/plot_spectra.py -s mutspec192.tsv -o ${label}_ms192.pdf -t

else
    echo "ArgumentError: nspecies"
    exit 1
fi

cd -
