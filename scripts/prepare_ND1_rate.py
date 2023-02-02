import re

from Bio import SeqIO
import pandas as pd

rates = pd.read_csv("data/MUS/mus_rates/ND1.rate", sep='\t', comment="#")

non_letters_positions = set()  # 1-BASED
for rec in SeqIO.parse("data/exposure/mus_nd1/raw_aln.fna", format="fasta"):
    seq = str(rec.seq)
    for i in range(0, len(seq), 3):
        codon = seq[i: i+3]
        if not re.match("[A-Za-z]{3}", codon):
            non_letters_positions.add(i+1)
            non_letters_positions.add(i+2)
            non_letters_positions.add(i+3)

print(sorted(non_letters_positions))

df = rates[~rates.Site.isin(non_letters_positions)]
df.to_csv("data/exposure/mus_nd1/ND1.rate", index=False, sep="\t")