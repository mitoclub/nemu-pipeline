#!/bin/bash

PATH_TO_NEMU=/home/kpotoh/nemu-pipeline/pipeline/Nemu-core/main.nf
PATH_TO_CONFIG=/home/kpotoh/nemu-pipeline/data/alisim/nemu_mam.config

INDIR=generations_mam
OUTDIR=spectra_reconstructed_mam
# treefile=/home/kpotoh/nemu-pipeline/data/alisim/external/mam_with_outgrp.nwk

MAX_NJOBS=800
SLEEP_TIME=900 # secs

COUNTER=1

for aln in $INDIR/*.fa.leaves; do
    let COUNTER++
    bn=`basename $aln .fa.leaves`
    workdir=$OUTDIR/$bn

    if [ -f $treefile ] && [ -f $aln ]; then
        # echo $bn
        mkdir -p $workdir
        cd $workdir
        # -resume
        nextflow -q -bg -c $PATH_TO_CONFIG run $PATH_TO_NEMU -with-trace \
            --sequence ../../$aln --outdir .
        cd -

        if [ `expr $COUNTER % $MAX_NJOBS` -eq 0 ]; then
            echo sleep $SLEEP_TIME secs
            sleep $SLEEP_TIME
        fi
    fi

done
