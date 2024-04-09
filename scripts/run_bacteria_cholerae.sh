#!/bin/bash

PATH_TO_NEMU=/home/kpotoh/nemu-pipeline/pipeline/nemu-core.nf
PATH_TO_CONFIG=/home/kpotoh/nemu-pipeline/data/cholerae/nemu_cholerae.config
PATH_TO_INPUT=/home/kpotoh/nemu-pipeline/data/cholerae/input
PATH_TO_OUTPUT=/home/kpotoh/nemu-pipeline/data/cholerae/output

if [ ! -f $PATH_TO_NEMU ]
then 
    echo "$PATH_TO_NEMU doesn't exist"
    exit 1
fi
if [ ! -f $PATH_TO_CONFIG ]
then 
    echo "$PATH_TO_CONFIG doesn't exist"
    exit 1
fi
if [ ! -e $PATH_TO_OUTPUT ]
then 
    echo "$PATH_TO_OUTPUT doesn't exist"
    exit 1
fi
if [ ! -e $PATH_TO_INPUT ]
then 
    echo "$PATH_TO_INPUT doesn't exist"
    exit 1
fi


MAX_NJOBS=64
SLEEP_TIME=900 # secs
RUNLIMIT=10
COUNTER=1

for species in `ls $PATH_TO_INPUT`; do
    if [ $COUNTER -gt $RUNLIMIT ]; then
        exit 0
    fi

    input_fasta=$PATH_TO_INPUT/$species
    workdir=$PATH_TO_OUTPUT/$(basename $species .fasta)

    if [ -e $workdir ]; then 
        echo "'$workdir' already exist: pass"
        continue
    fi
    let COUNTER++

    echo -e "workdir: $workdir"
    mkdir -p $workdir
    cd $workdir
    cp $input_fasta .

    aln=$species

    # -resume
    # -q -bg
    nextflow -bg -q -c $PATH_TO_CONFIG run $PATH_TO_NEMU -with-trace --sequence $aln --outdir .
    cd - >/dev/null
    sleep 0.5

    if [ `expr $COUNTER % $MAX_NJOBS` -eq 0 ]; then
        echo sleep $SLEEP_TIME secs
        sleep $SLEEP_TIME
    fi
done
echo

echo "Processed $COUNTER fastas"  # 17748

# rm -rf $PATH_TO_OUTPUT/*/*/*/work
