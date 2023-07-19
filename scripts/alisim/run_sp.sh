#!/bin/bash

cd /home/kpotoh/nemu-pipeline/data/alisim

PATH_TO_NEMU=/home/kpotoh/nemu-pipeline/pipeline/Nemu-core/main.nf
PATH_TO_CONFIG=/home/kpotoh/nemu-pipeline/data/alisim/nemu_sp.config

INDIR=generations_sp
OUTDIR=spectra_reconstructed_sp

MAX_NJOBS=40
SLEEP_TIME=10 # secs; single sp100 processed within 170sec

COUNTER=1

for treefile in $INDIR/*.nwk; do
    let COUNTER++
    aln=${treefile%.nwk}.fa.leaves
    bn=`basename $treefile .nwk`
    workdir=$OUTDIR/$bn

    if [ -f $treefile ] && [ -f $aln ]; then
        # echo $bn
        mkdir -p $workdir
        cd $workdir
        # -resume
        nextflow -q -bg -c $PATH_TO_CONFIG run $PATH_TO_NEMU -with-trace -resume  \
            --sequence ../../$aln --outdir . --treefile ../../$treefile --verbose false
        cd - >/dev/null

        if [ `expr $COUNTER % $MAX_NJOBS` -eq 0 ]; then
            echo sleep $SLEEP_TIME secs
            sleep $SLEEP_TIME
        fi
    fi

done
