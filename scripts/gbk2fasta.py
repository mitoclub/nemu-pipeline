#!/bin/bash
# USAGE: gbk2fasta.py [gbk_files...] out_fasta

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

fout = open(sys.argv[-1], "w")
for i, file in enumerate(sys.argv[1:-1]):
    rec: SeqRecord = None
    for rec in SeqIO.parse(file, "genbank"):
        try:
            if rec.seq is not None:
                SeqIO.write([rec], fout, "fasta")
        except:
            print("pass {}".format(rec.name))
fout.close()
