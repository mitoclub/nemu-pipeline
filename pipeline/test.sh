#!/bin/bash

if [ -f nemu-core.nf ] && [ -f test-core.config ]; then
    seq=../example_data/CO1__Stenella_longirostris.fna
    nextflow -c test-core.config run nemu-core.nf --sequence $seq --outdir test_output1
else
    echo you must be in the pipeline directory to access nemu-core.nf and test-core.config files
    exit 1
fi

# if [ -f nemu.nf ] && [ -f test.config ]; then
#     seq=../example_data/CYTB__Stenella_longirostris.faa
#     # nextflow -c test.config run nemu.nf --sequence $seq --outdir test_output2
#     nextflow -c test.config run nemu.nf --sequence $seq --outdir test_output2
# else
#     echo you must be in the pipeline directory to access nemu.nf and test.config files
#     exit 1
# fi

echo DONE
