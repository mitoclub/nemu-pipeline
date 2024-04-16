#!/usr/bin/env python3
# USAGE: script_name.py query_msa.fasta [mutations.tsv]

import sys
from functools import reduce, partial
from multiprocessing import Pool

import pandas as pd
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.AlignInfo import SummaryInfo
from pymutspec.annotation import CodonAnnotation

import warnings
warnings.filterwarnings('ignore')

coda = CodonAnnotation(gencode=2)


def read_ingroup_msa(path):
    msa = AlignIO.read(path, 'fasta')
    for i,x in enumerate(msa):
        if x.name == 'OUTGRP':
            break
    assert msa[i].name == 'OUTGRP'
    del msa[i]
    return msa


def read_ingroup_msa_and_outgrp(path):
    msa = AlignIO.read(path, 'fasta')
    for i,x in enumerate(msa):
        if x.name == 'OUTGRP':
            break
    assert msa[i].name == 'OUTGRP'
    outgrp = str(msa[i].seq)
    del msa[i]
    return msa, outgrp


def collect_snp_from_consensus(rec: SeqRecord, consensus: str):
    seq = str(rec.seq).upper()
    # print(seq)
    # for i, (n1, n2) in enumerate(zip(seq, consensus)):
        # if n1 != n2:
        #     print('site', i, n1, n2, seq[i-2:i+3], consensus[i-2:i+3])
    muts = coda.extract_mutations_simple(consensus, seq)
    return muts


def milti_collect(x, consensus):
    i, rec = x
    m = collect_snp_from_consensus(rec, consensus)
    if len(m):
        m = m.assign(seq_id=i)
    return m


# p = Pool(64)
# _ = p.map(partial(milti_collect, consensus=consensus), enumerate(msa))
# x = pd.concat(_, ignore_index=True)


def main():
    if sys.argv[1] == '-h' or sys.argv[1] == '-help' or sys.argv[1] == '--help':
        print(f"USAGE: {sys.argv[0]} query_msa.fasta [mutations.tsv]", file=sys.stderr)
        exit(0)
    elif sys.argv[1] == '-h' or sys.argv[1] == '-help' or sys.argv[1] == '--help':
        print(f"USAGE: {sys.argv[0]} query_msa.fasta [mutations.tsv]", file=sys.stderr)
        raise RuntimeError('Please specify input fasta file with nucleotide alignment')
    elif len(sys.argv) == 2:
        path_fasta = sys.argv[1]
        outpath = None
    elif len(sys.argv) >= 3:
        path_fasta = sys.argv[1]
        outpath = sys.argv[2]

    msa = read_ingroup_msa(path_fasta)
    si = SummaryInfo(msa)
    consensus = si.dumb_consensus(0.5, 'N').upper()
    # print(consensus)
    
    data = []
    for i, rec in enumerate(msa):
        # print('seq', i)
        m = collect_snp_from_consensus(rec, consensus)
        # print(m)
        if len(m):
            data.append(m.assign(seq_id=i))
    if len(data):
        mut = pd.concat(data, ignore_index=True)
        if outpath:
            mut.to_csv(outpath, index=False, sep='\t')
        else:
            print(mut)
    else:
        print("No mutations found", file=sys.stderr)

if __name__ == '__main__':
    main()

    # print("TEST")
    # print(coda.extract_mutations_simple("ATCGAC", "ACCGTC"))