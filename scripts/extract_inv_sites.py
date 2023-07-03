#!/usr/bin/env python3

from collections import defaultdict
import click
import pandas as pd
from Bio import SeqIO

@click.command("extrctor", help="")
@click.argument("infile", type=click.Path(True, True), required=True)
@click.argument("outfile", type=click.Path(False), required=True)
# @click.option("--aln", "aln_dir", required=True, type=click.Path(True), help="Path to directory with gene alignment files. Used for gap filling")
def main(infile, outfile):
    site_rates = defaultdict(lambda: defaultdict(int))
    fasta = SeqIO.parse(infile, "fasta")
    for rec in fasta:
        seq = str(rec.seq)
        for site, nuc in enumerate(seq, 1):
            site_rates[site][nuc] += 1

    d = pd.DataFrame(site_rates).T
    nseqs = d.max().max()
    inv = (d.max(axis=1) != nseqs).astype(int).rename("Cat")
    inv.index.name = "Site"
    inv.to_csv(outfile, sep="\t")


if __name__ == "__main__":
    main()
