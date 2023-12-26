#!/bin/bash

cd /home/kpotoh/nemu-pipeline/data/bacteria

PATH_TO_NEMU=/home/kpotoh/nemu-pipeline/pipeline/nemu-core.nf
PATH_TO_CONFIG=/home/kpotoh/nemu-pipeline/data/bacteria/nemu_bacteria.config
PATH_TO_INPUT=/scratch/genkvg/bacteria/askudnov
PATH_TO_OUTPUT=/scratch/bacteria/spectra

# _output1: 1786
# _output2: 6377
# _output3: 5045
CLADES=("_output1" "_output2" "_output3")

MAX_NJOBS=60
SLEEP_TIME=120 # secs
RUNLIMIT=12000
COUNTER=1

for clade in ${CLADES[@]}; do
    echo $clade
    for species in `ls $PATH_TO_INPUT/$clade/ | grep _`; do
        INDIR=$PATH_TO_INPUT/$clade/$species
        OUTDIR=$PATH_TO_OUTPUT/$clade/$species

        for aln in $INDIR/*.fasta; do
            if [ $COUNTER -gt $RUNLIMIT ]; then
                exit 0
            fi

            workdir=$OUTDIR/$(basename $aln .fasta)
            if [ -e $workdir ]; then 
                echo "'$workdir' already exist: pass"
                continue
            fi
            let COUNTER++

            aln_abs=`realpath $aln`

            echo -e "workdir: $workdir"
            mkdir -p $workdir
            cd $workdir

            # -resume
            # -q -bg
            nextflow -bg -q -c $PATH_TO_CONFIG run $PATH_TO_NEMU -with-trace --sequence $aln_abs --outdir .
            cd - >/dev/null
            sleep 0.5

            if [ `expr $COUNTER % $MAX_NJOBS` -eq 0 ]; then
                echo sleep $SLEEP_TIME secs
                sleep $SLEEP_TIME
            fi
        done
        echo
    done
    echo -e '\n'
done

echo "Processed $COUNTER fastas"  # 17748

# rm -rf $PATH_TO_OUTPUT/*/*/*/work
