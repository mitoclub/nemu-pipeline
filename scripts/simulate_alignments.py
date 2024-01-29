#!/usr/bin/env python3
'''
Simulation of multiple nucleotide alignments based on 
nucleotide alignment, tree file, 12-comp spectrum (rates) and site rates
'''

from collections import defaultdict
import random
from typing import Dict

import click
import pyvolve
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pymutspec.io import read_rates

DEFAULT_REPLICS = 10
DEFAULT_GENCODE = 2
DEFAULT_RATECATCUTOFF = 1
EPSILON = 1e-10


def read_rates(path_to_mutspec, eps=EPSILON):
    ms = pd.read_csv(path_to_mutspec, sep="\t")
    ms["Mut"] = ms["Mut"].str.replace(">", "")
    ms["MutSpec"] = ms["MutSpec"] + eps
    rates = ms.set_index("Mut")["MutSpec"].to_dict()
    return rates


def get_root_seq(path_to_fasta):
    root_seq = None
    for rec in SeqIO.parse(path_to_fasta, format="fasta"):
        if rec.name.startswith("RN"):
            root_seq = str(rec.seq)
            break
    return root_seq


def get_root_cdn(path_to_fasta, cdn_pos: int, r=0.015, verbose=False):
    '''select random sequence from fasta file. Condition: get 1 random seq from 100
    
    cdn_pos - 1-based'''
    i = 0
    for rec in SeqIO.parse(path_to_fasta, format="fasta"):
        i += 1
        if rec.name.startswith("RN"):
            nuc_start = (cdn_pos - 1) * 3 
            nuc_end = nuc_start + 3
            cdn = str(rec.seq[nuc_start: nuc_end])
            if r > random.random():
                if verbose:
                    print(f'seq id: {i}')
                break
    return cdn


def codonify(seq: str):
    codons = []
    for i in range(0, len(seq), 3):
        codons.append(seq[i: i + 3])
    return codons


def seq_masking(seq, codons):
    """
    Arguments
    ---------
    seq: str
        root seq
    codons: array
        codon mask, 1-based
    """
    if len(seq) // 3 < max(codons):
        raise ValueError("codon mask don't fit to sequence")

    seq_cdn = codonify(seq)
    masked_codons = [seq_cdn[i - 1] for i in codons]
    return "".join(masked_codons)


def seq_unmasking(seq, codons, consensus):
    if len(consensus) // 3 < max(codons):
        raise ValueError("codon mask don't fit to sequence")

    seq_cdn = codonify(seq)
    consensus = codonify(consensus)
    codons_set = set(codons)
    unmasked = ''
    var_cdn_id = 0
    for i, cdn in enumerate(consensus, 1):
        if i in codons_set:
            cdn_var = seq_cdn[var_cdn_id]
            unmasked += cdn_var
            var_cdn_id += 1
        else:
            unmasked += cdn
    return unmasked


def msa_unmasking(consensus, codons, aln_dct: Dict[str, str]) -> Dict[str, str]:
    """
    Arguments
    ---------
    consensus: str
        root or some seq to gen conservative positions from
    codons: array
        codon mask, 1-based
    aln: List[str]
        simulated sequences without masked codons
    """
    if len(consensus) // 3 < max(codons):
        raise ValueError("codon mask don't fit to sequence")

    codons_set = set(codons)
    msa_rebuilt = dict()
    for node, seq in aln_dct.items():
        msa_rebuilt[node] = seq_unmasking(seq, codons_set, consensus)
    return msa_rebuilt


def write_seqs(seqdict: Dict[str, str], seqfile, seqfmt="fasta-2line"):
    alignment = []
    for entry in seqdict:
        seq_entry = Seq(seqdict[entry])
        seq_object = SeqRecord(seq_entry, id = entry, description = "", annotations={"molecule_type": "DNA"}) 
        alignment.append(seq_object)     
    SeqIO.write(alignment, seqfile, seqfmt)


def get_variable_codon_positions(ratesfile, ratecat_cutoff=1):
    mask = (read_rates(ratesfile) > ratecat_cutoff).astype(np.int8)
    mask = mask[:len(mask) - len(mask)%3]
    all_codons = np.reshape(mask, (-1, 3))
    # indexes of codons that changed, 1-based
    changed_codons_positions = np.where(all_codons.sum(axis=1) > 0)[0] + 1
    return changed_codons_positions


def _get_codon_freqs_looong(filepath, columns, gencode):
    reader = pyvolve.ReadFrequencies("codon", file=filepath, gencode=gencode, columns=columns)
    codonfreqs = reader.compute_frequencies(type="codon")
    codonfreqs += EPSILON
    codonfreqs /= np.sum(codonfreqs)
    return codonfreqs


def get_codon_freqs(filepath, codons: list=None, gencode=2):
    '''
    Arguments
    ------
    - filepath, str - path to nucleotide alignment
    - codons, List[int] - list of codon positions (1-based) to count. Not the same as nuc position
    - gencode, int

    Return
    ------
    - freqs, dict[int, list[float]] - dict of lists with codon freqs for codon positions
    '''
    counts = defaultdict(lambda: defaultdict(int)) # {cdn_pos: {cdn: count}}
    msa = SeqIO.parse(filepath, 'fasta')
    for rec in msa:
        seq_cdn = codonify(str(rec.seq))
        if codons is None:
            codons = range(1, len(seq_cdn) + 1)
        for cdn_pos in codons: # 1-based
            cdn = seq_cdn[cdn_pos - 1]
            counts[cdn_pos][cdn] += 1
    
    g = pyvolve.Genetics(gencode)
    freqs = {k: [v[c]/sum(v.values())+EPSILON for c in g.codons] for k,v in counts.items()}
    freqs = {k: [f/sum(v) for f in v] for k,v in freqs.items()}
    return freqs


@click.command("MutSel simulation", help="")
@click.option("-a", "--alignment", "path_to_mulal", type=click.Path(True), help="")
@click.option("-t", "--tree", "path_to_tree", type=click.Path(True), help="")
@click.option("-s", "--spectra", "path_to_mutspec", type=click.Path(True), help="")
@click.option("--rates", type=click.Path(True), default=None, help="path to rates from iqtree that will be used for positions masking")
@click.option("-o", "--outseqs", type=click.Path(writable=True), help="Output sequences alignment (fasta)")
@click.option("--outcount", type=click.Path(writable=True), default=None, help="")
@click.option("-r", "--replics", "nreplics", default=DEFAULT_REPLICS, show_default=True, type=int, help="")
@click.option("-w", "--write_anc", is_flag=True, help="")
@click.option("-c", "--gencode", default=DEFAULT_GENCODE, show_default=True, help="")
@click.option("--ratecat_cutoff", default=DEFAULT_RATECATCUTOFF, show_default=True, help="")
# @click.option("-l", "--scale_tree", default=1., show_default=True, help="")
def main(path_to_mulal, path_to_tree, path_to_mutspec, outseqs, outcount, nreplics, write_anc, gencode, ratecat_cutoff, rates):
    tree = pyvolve.read_tree(file=path_to_tree)  # scale_tree=scale_tree_factor
    custom_mutation_asym = read_rates(path_to_mutspec)
    codons = get_variable_codon_positions(rates, ratecat_cutoff) if rates else None
    print(f'rates: {custom_mutation_asym}')
    print(f'#variable codons: {len(codons)}')
    freqs = get_codon_freqs(path_to_mulal, codons, 2)
    consensus = get_root_seq(path_to_mulal) # first seq from alignment

    partitions = []
    for pos in codons:
        mutsel_params = {
            "state_freqs": freqs[pos], 
            "mu": custom_mutation_asym
        }
        model = pyvolve.Model("mutsel", mutsel_params, gencode=gencode)
        part = pyvolve.Partition(models=model, size=1)
        partitions.append(part)

    evolver = pyvolve.Evolver(partitions=partitions, tree=tree, gencode=gencode)

    for i in range(nreplics):
        print("Generating replica {}".format(i))
        seqfile = outseqs.replace(".fasta", f"_sample-{i:04}.fasta")
        cur_outcount = outcount + f"_sample-{i:04}"
        if codons is None:
            evolver(
                seqfile=seqfile, 
                countfile=cur_outcount,
                ratefile=None, infofile=None,
                write_anc=write_anc,
            )
        else:
            evolver(
                seqfile=None,
                countfile=cur_outcount,
                ratefile=None, infofile=None,
                write_anc=write_anc,
            )
            msa_dct = evolver.get_sequences(True)
            unmasked_aln = msa_unmasking(consensus, codons, msa_dct)
            print(unmasked_aln)


            # write_seqs(unmasked_aln, seqfile)
    print('Done')


if __name__ == "__main__":
    gene = 'mouse_cytb'
    folder = f'./data/selection_search/{gene}'

    msa = f'{folder}/pyvolve/mulal.fasta.clean'
    tree = f'{folder}/pyvolve/tree.nwk.ingroup'
    sp = f'{folder}/ms/ms12syn_internal_{gene}.tsv'
    sr = f'./data/selection_search/rates/{gene}.rate'
    out = f'test_simulation/{gene}.fasta'
    outcount='test_simulation/outcount'
    # main(f"-a {msa} -t {tree} -s {sp} -o {out} --rates {sr} -w -r 2 -c 2 --outcount {outcount}".split())
    
    import json
    a = open('consensus.txt')
    consensus = json.load(a)
    a.close()
    a = open('codons.txt')
    codons = json.load(a)
    a.close()
    a = open('aln_dct.txt')
    aln_dct = json.load(a)
    a.close()
    
    freqs = get_codon_freqs(msa, codons, 2)
    x = msa_unmasking(consensus, codons, aln_dct)

    codons = [1,3,5]
    seqs =     ['ACAAACATACGAAAA']
    consensus = 'XXXXXXXXXXXXXXX'
    print(codons)

    # seqs = [str(x.seq) for x in SeqIO.parse(msa, 'fasta')]
    # seqs_masked = [seq_masking(x, codons) for x in seqs]
    # seqs_unmasked = [seq_unmasking(x, codons, consensus) for x in seqs_masked]
    # print(codonify(seqs[0]))
    # print(codonify(seqs_masked[0]))
    # print(codonify(seqs_unmasked[0]))


    # for i, seq in enumerate(seqs):
    #     seq_unmasked = codonify(seqs_unmasked[i])
    #     seq = codonify(seq)
    #     for j, cdn in enumerate(seq):
    #         cdn_unmasked = seq_unmasked[j]
    #         if cdn != cdn_unmasked:
    #             print(f'{i}, {j}, {cdn}, {cdn_unmasked}')
        # break
    
    # assert seqs == seqs_unmasked

    # main()
    # codons = get_variable_codon_positions(path_to_rates, 1)
    # print(f'#variable codons: {len(codons)}')
    # freqs = get_codon_freqs(path_to_mulal, codons, 2)
    # print(freqs)