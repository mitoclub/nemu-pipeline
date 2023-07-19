#!/bin/bash
# USAGE: script.sh NWK OUTDIR

PATH_TO_NEMU=/home/kpotoh/nemu-pipeline/pipeline/Nemu-core/main.nf
# PATH_TO_CONFIG=/home/kpotoh/nemu-pipeline/data/alisim/nemu_sp.config
PATH_TO_CONFIG=/home/kpotoh/nemu-pipeline/data/alisim/nemu_mam.config


if [ ! -e $2 ]; then
    echo "outdir '$2' doesn't exist"
    exit 1;
fi

if [ ! -f $1 ]; then
    echo "input tree '$1' doesn't exist"
    exit 1;
fi

treefile=$1
OUTDIR=$2

aln=${treefile%.nwk}.fa.leaves
bn=`basename $treefile .nwk`
workdir=$OUTDIR/$bn
if [ -f $treefile ] && [ -f $aln ]; then
    mkdir -p $workdir
    cd $workdir
    nextflow -c $PATH_TO_CONFIG run $PATH_TO_NEMU -with-trace \
        --sequence ../../$aln --outdir . --treefile ../../$treefile
    cd -
else
    echo "Input files ($treefile, $aln) doesn't exist"
fi
