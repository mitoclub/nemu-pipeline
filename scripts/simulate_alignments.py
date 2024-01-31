#!/usr/bin/env python3
'''
Simulation of multiple nucleotide alignments based on 
nucleotide alignment, tree file, 12-comp spectrum (rates) and site rates
'''

import copy
import sys
import random
from collections import defaultdict
from typing import Dict

import click
import pyvolve
import tqdm
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pymutspec.io import read_rates
from traitlets import default

EPSILON = 1e-10


def read_spectrum(path_to_mutspec, eps=EPSILON):
    ms = pd.read_csv(path_to_mutspec, sep="\t")
    ms["Mut"] = ms["Mut"].str.replace(">", "")
    ms["MutSpec"] = ms["MutSpec"] + eps
    rates = ms.set_index("Mut")["MutSpec"].to_dict()
    return rates


def codonify(seq):
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


def msa_unmasking(consensus, codons, msa_dct: Dict[str, str]) -> Dict[str, str]:
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
    for node, seq in msa_dct.items():
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
    mask = mask[:len(mask) - len(mask) % 3]
    all_codons = np.reshape(mask, (-1, 3))
    # indexes of codons that changed, 1-based
    changed_codons_positions = np.where(all_codons.sum(axis=1) > 0)[0] + 1
    return changed_codons_positions


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


def get_consensus(path_to_msa, gencode):
    freqs = get_codon_freqs(path_to_msa, None, gencode)
    g = pyvolve.Genetics(gencode)
    n = max(freqs.keys())
    consensus = []
    for i in range(1, n + 1):
        site_freqs = freqs[i]
        pos_max = np.argmax(site_freqs)
        cdn = g.codons[pos_max]
        consensus.append(cdn)
    return ''.join(consensus)


def scale_tree(t: pyvolve.newick.Node, factor=1):
    if isinstance(t.branch_length, (int, float)):
        t.branch_length *= factor
    if len(t.children) == 2:
        c1, c2 = t.children
        scale_tree(c1, factor)
        scale_tree(c2, factor)
    elif len(t.children) == 0:
        return
    else:
        raise Exception('Not a binary tree!')


def read_codon_specific_trees(path_to_tree, path_to_rec_mut, msa_len: int, par_factor):
    mut = pd.read_csv(path_to_rec_mut, sep='\t')
    print(f'nmut (all): {len(mut)}', file=sys.stderr)
    mut = mut[mut.Label > 0]  # only syn
    print(f'nmut (syn): {len(mut)}', file=sys.stderr)
    assert mut.PosInGene.max() <= msa_len
    codons = range(1, msa_len // 3 + 1) # 1-based
    mean_mut_rate = len(mut) / msa_len
    tree_base = pyvolve.read_tree(file=path_to_tree)
    trees = dict()
    for pos in tqdm.tqdm(codons, 'Scale trees'):
        nuc_start = (pos - 1) * 3 + 1 # 1-based
        nuc_end = nuc_start + 3
        muts_in_codon = mut[mut.PosInGene.between(nuc_start, nuc_end, 'left')]
        cdn_site_mut_rate = len(muts_in_codon) / 3
        scale_tree_factor = cdn_site_mut_rate / mean_mut_rate * par_factor + EPSILON
        tree = copy.deepcopy(tree_base)
        scale_tree(tree, scale_tree_factor)
        trees[pos] = tree
    return trees


def read_parsimony_informative_sites_num(path_to_log):
    with open(path_to_log, 'r') as fin:
        for line in fin:
            if 'Number of parsimony informative sites:' in line:
                n = line.strip().split()[-1]
                try:
                    n = int(n)
                except Exception as e:
                    raise e
                break
    return n


def run_simulation(trees, freqs, mutation_rates, codons, gencode=2):
    msa_dct = defaultdict(str)
    for pos in tqdm.tqdm(codons, 'Simulate codons'):
        mutsel_params = {"state_freqs": freqs[pos], "mu": mutation_rates}
        model = pyvolve.Model("mutsel", mutsel_params, gencode=gencode)
        part = pyvolve.Partition(models=model, size=1)
        evolver = pyvolve.Evolver(partitions=part, tree=trees[pos], gencode=gencode)
        
        # run simulation
        evolver(
            seqfile=None, countfile=None,
            ratefile=None, infofile=None, write_anc=True,
        )
        cdn = evolver.get_sequences(True)
        for node, x in cdn.items():
            msa_dct[node] += x
    
    return msa_dct


@click.command("Simulation", help="")
@click.option("-a", "--alignment", "path_to_msa", type=click.Path(True), help="")
@click.option("-t", "--tree", "path_to_tree", type=click.Path(True), help="")
@click.option("-s", "--spectra", "path_to_mutspec", type=click.Path(True), help="")
@click.option("-o", "--outseqs", type=click.Path(writable=True), help="Output sequences alignment (fasta)")
@click.option("--rates", type=click.Path(True), default=None, help="path to rates from iqtree that will be used for positions masking")
@click.option("--rec", "path_to_rec", type=click.Path(True), help='path to table with reconstructed (observed) mutations')
@click.option("--iqtree-log", type=click.Path(True), default=None, help='path to iqtree_anc_report.log')
@click.option("-r", "--replics", default=10, show_default=True, type=int, help="")
@click.option("-c", "--gencode", default=2, show_default=True, help="")
@click.option("--ratecat_cutoff", default=1, show_default=True, help="")
def main(path_to_msa, path_to_tree, path_to_mutspec, rates, path_to_rec, iqtree_log, outseqs, replics, gencode, ratecat_cutoff):
    custom_mutation_asym = read_spectrum(path_to_mutspec)
    codons = get_variable_codon_positions(rates, ratecat_cutoff) if rates else None
    consensus = get_consensus(path_to_msa, gencode)
    used_codons_share = len(codons) * 3 / len(consensus) * 100
    print(f'rates: {custom_mutation_asym}\n#variable codons: {len(codons)}\nused codons share: {used_codons_share:.2f}%', file=sys.stderr)
    freqs = get_codon_freqs(path_to_msa, codons, 2)
    n_par_sites = read_parsimony_informative_sites_num(iqtree_log) if iqtree_log else len(codons)
    par_factor = n_par_sites // 3
    trees = read_codon_specific_trees(path_to_tree, path_to_rec, len(consensus), par_factor)

    for i in range(replics):
        print(f"Generating replica {i}", file=sys.stderr)
        seqfile = outseqs.replace(".fasta", f"_sample-{i:04}.fasta")
        msa_dct = run_simulation(trees, freqs, custom_mutation_asym, codons, gencode)
        unmasked_aln = msa_unmasking(consensus, codons, msa_dct)
        write_seqs(unmasked_aln, seqfile)
    print('Done', file=sys.stderr)


'''Надо будет ещё пересмотреть последующий анализ мутаций. 
Мы сравниваем конкретные количества мутаций без нормировки, но чтобы так делать, нужно получить хотя бы похожие количества мутаций в реконструкции и симуляции. 
Нам ещё не удалось приблизиться к такому,'''


if __name__ == "__main__":
    main()
    # gene = 'mouse_cytb'
    # folder = f'./data/selection_search/{gene}'

    # msa = f'{folder}/pyvolve/mulal.fasta.clean'
    # tree = f'{folder}/pyvolve/tree.nwk.ingroup'
    # sp = f'{folder}/ms/ms12syn_internal_{gene}.tsv'
    # sr = f'./data/selection_search/rates/{gene}.rate'
    # rec = f'{folder}/tables/observed_mutations_iqtree.tsv'
    # llog = f'{folder}/logs/iqtree_anc_report.log'
    # out = f'test_simulation/{gene}.fasta'
    # main(f"-a {msa} -t {tree} -s {sp} -o {out} --rates {sr} --rec {rec} --iqtree-log {llog} -r 2 -c 2".split())


# collect_mutations.py --tree data/selection_search/mouse_cytb/pyvolve/tree.nwk \
# --states test_simulation/mouse_cytb_sample-0000.fasta --states-fmt fasta \
# --gencode 2 --syn --no-mutspec --outdir test_simulation/out0 --force
