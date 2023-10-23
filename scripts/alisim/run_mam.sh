#!/bin/bash

# PATH_TO_NEMU=/home/kpotoh/nemu-pipeline/data/alisim/nemu.nf
# PATH_TO_CONFIG=/home/kpotoh/nemu-pipeline/data/alisim/nemu_mam.config
PATH_TO_NEMU=/home/kpotoh/nemu-pipeline/pipeline/Nemu-core/nemu.nf
PATH_TO_CONFIG=/home/kpotoh/nemu-pipeline/data/alisim/nemu_mam.config

INDIR=generations_mam2
OUTDIR=spectra_reconstructed_mam2
# treefile=/home/kpotoh/nemu-pipeline/data/alisim/external/mam_with_outgrp.nwk

MAX_NJOBS=52
SLEEP_TIME=4000 # secs
WAIT_TIME=20

COUNTER=1

# for aln in $INDIR/12.12*.fa.leaves; do
for aln in $INDIR/*.fa.leaves; do
    let COUNTER++
    bn=`basename $aln .fa.leaves`
    workdir=$OUTDIR/$bn

    if [ -f $treefile ] && [ -f $aln ]; then
        # echo $bn
        mkdir -p $workdir
        cd $workdir
        # -q -bg   run -resume 
        nextflow -q -bg -c $PATH_TO_CONFIG run -resume $PATH_TO_NEMU -with-trace \
            --sequence ../../$aln --outdir .
        cd -
        sleep $WAIT_TIME

        if [ `expr $COUNTER % $MAX_NJOBS` -eq 0 ]; then
            echo sleep $SLEEP_TIME secs
            sleep $SLEEP_TIME
        fi
    fi

done
