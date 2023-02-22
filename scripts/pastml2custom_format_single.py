#!/usr/bin/env python3

import os
import re
import glob
import sys

import click
import numpy as np
import pandas as pd


def read_marginal_probabilities(filepath, model="F81"):
    nucleotides = list("ACGT")
    character = re.search(f"marginal_probabilities\.character_(.+)\.model_{model}\.tab", filepath)
    states = pd.read_csv(filepath, sep="\t")
    states.rename(columns={"node": "Node"}, inplace=True)
    states["State"] = [states.columns[i] for i in (np.argmax(states.iloc[:, 1:].values, axis=1) + 1)]
    for nucl in nucleotides:
        if nucl in states.columns:
            states["p_" + nucl] = states[nucl].astype(np.float16)
            states.drop(nucl, axis=1, inplace=True)
        else:
            states["p_" + nucl] = np.float16(0.0)
    states["Site"] = int(character.groups()[0])
    states["Part"] = 1
    states = states[['Node', 'Part', 'Site', 'State', 'p_A', 'p_C', 'p_G', 'p_T']]
    return states


def read_data(indir, model: str):
    """
    Arguments
    ------
    indir: str
        path to dir, containing files with states for each separate character
    """
    gene = []
    for path in glob.glob(os.path.join(indir, "*")):
        if not "marginal_probabilities" in path:
            continue
        df = read_marginal_probabilities(path, model)
        gene.append(df)
    gene_df = pd.concat(gene)
    return gene_df


@click.command("formatter", help="Reformat pastml output to usual states table")
@click.option("-i", "--indir", type=click.Path(True, False), required=True, help="Path to the directory that contains pastml output files")
@click.option("-o", "--outpath", required=True, type=click.Path(writable=True), help="Path to output states file (tsv)")
@click.option("-m", "--model", type=str, default="F81", show_default=True, help="Model used to calculate ancestral states")
def main(indir, outpath, model):
    gene_states = read_data(indir, model)
    print("Sorting...", file=sys.stderr)
    gene_df = gene_states.sort_values(["Node", "Part", "Site"])
    print("Writing...", file=sys.stderr)
    if gene_df.Node.value_counts().unique().shape[0] > 1:
        print("WARNING: some positions absent in some sequences!!!")
        print(gene_df.Node.value_counts().sort_values())
    gene_df.to_csv(outpath, index=None, sep="\t")


if __name__ == "__main__":
    main()

