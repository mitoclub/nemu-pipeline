import os
import sys
import random
from collections import Counter

from scipy.stats import entropy
from Bio import SeqIO
from ete3 import PhyloTree

MAMMALS_TREE  = "/home/kpotoh/nemu-pipeline/data/alisim/external/mam_with_outgrp.nwk"
ROOT_SEQ      = "/home/kpotoh/nemu-pipeline/data/alisim/external/human_cytb.fasta"
ROOT_SEQ_NAME = "human_cytb"
LENGTH = 1140
SEED = random.randint(0, 100000)

debug_mode = True
OUTDIR = "generations_mam"
N = 50

nucls = list('ACGT')
base = 4
human_cytb_seq = str(next(SeqIO.parse(ROOT_SEQ, 'fasta')).seq)
human_cytb_freqs = dict(Counter(human_cytb_seq))
human_cytb_freqs = [human_cytb_freqs[x]/len(human_cytb_seq) for x in nucls]
human_cytb_entropy = entropy(human_cytb_freqs, base=base)


cmd = "iqtree2 --alisim {} -m {}{}{}+G6+I{} -t {} --seed {} --write-all {} -af fasta {} -redo"


if not os.path.exists(OUTDIR):
    os.mkdir(OUTDIR)
    print(f"{OUTDIR} created", file=sys.stderr)


def beautify_array(lst):
    lst = [str(round(x, 2)) for x in lst]
    return "{" + "/".join(lst) + "}"


def generate_f(as_natural=False) -> list:
    while True:
        a = [random.random() for _ in range(4)]
        s = sum(a)
        a = [x / s for x in a]
        if as_natural:
            entr = entropy(a, base=base)
            if abs(human_cytb_entropy - entr) < 0.05:
                break
        else:
            break
    return a


def generate_rates(n, rev=True):
    if rev:
        r = [random.lognormvariate(-2, 1) for _ in range(n)]
    else:
        param_max = 0.97
        r = [random.normalvariate(0, param_max / 2) for _ in range(n)]
        r = [max(min(x, param_max), -param_max) for x in r]
    
    return r


def add_outgrp(tree: PhyloTree, factor=5, inplace=True):
    '''just change T1 node distance and rename it'''
    if not inplace:
        tree = tree.copy()
    leaves_dists = [x.dist for x in tree.iter_leaves()]
    d = sum(leaves_dists) / len(leaves_dists) * factor
    outgrp = tree.search_nodes(name="T1")[0]
    outgrp.dist = d
    outgrp.name = "OUTGRP"
    tree.set_outgroup("OUTGRP")
    tree.name = f"Node{tree.name}"
    for innode in tree.iter_descendants():
        if not innode.is_leaf() and len(innode.name) > 0:
            innode.name = f"Node{innode.name}"
    return tree


def generate_mdl_params(mdl: str):
    if mdl == "gtr":
        mdl_params = beautify_array(generate_rates(6, True))
    elif mdl == "12.12":
        mdl_params = beautify_array(generate_rates(11, False))
    elif mdl == "10.12":
        mdl_params = beautify_array(generate_rates(9, False))
    else:
        raise NotImplementedError
    return mdl_params


def main():
    log_file = open(f'alisim_generation_seed{SEED}.log', 'w')

    # mdl = "gtr"
    for mdl in ["gtr", "12.12"]:
        # for gene in ["rnd", "cytb"]:
        for gene in ["rnd",]:
            root_seq = f"--length {LENGTH}" if gene == "rnd" else f"--root-seq {ROOT_SEQ},{ROOT_SEQ_NAME}"

            # for raw_tree in [100, 1000]:
            for raw_tree in ["mam"]:
                if isinstance(raw_tree, int):
                    is_random_tree = True
                    tree = "RANDOM{bd{0.1/0.05}/" + str(raw_tree) + "}"
                    rlen_args = "-rlen 0 0.001 0.01"  # for intra-species trees with small branch lengths
                elif raw_tree == "mam":
                    is_random_tree = False
                    tree = MAMMALS_TREE
                    rlen_args = ""
                else:
                    raise NotImplementedError

                for i in range(N):
                    prefix = f"{OUTDIR}/{mdl}_{raw_tree}_{gene}_replica_{i}"
                    # mdl_params = "{0.5/0.6/0.9/0.2/0.1/0.4/0.7/0.8/0.3/0.15/0.65}"
                    mdl_params = generate_mdl_params(mdl)
                    f_custom = "+F" + beautify_array(generate_f(not is_random_tree)) if gene == "rnd" else ""
                    inv_rate = beautify_array([random.random() * 0.5])  # A proportion of invariable sites
                    cur_cmd  = cmd.format(prefix, mdl, mdl_params, f_custom, inv_rate, tree, SEED, root_seq, rlen_args)
                    log_file.write(cur_cmd + '\n')
                    if debug_mode:
                        print(cur_cmd)
                    else:
                        os.system(cur_cmd)

                    if is_random_tree and not debug_mode:
                        tree_without_outgrp = PhyloTree(prefix + ".full.treefile", format=1)
                        add_outgrp(tree_without_outgrp, random.randint(20, 30))
                        tree_without_outgrp.write(outfile=prefix + ".nwk", format=1)

                        for suffix in (".fa", ".treefile", ".full.treefile"):
                            os.remove(prefix + suffix)
                        
                        # generate again using trees with outgrp
                        pregenerated_tree = prefix + ".nwk"
                        cur_cmd = cmd.format(prefix, mdl, mdl_params, f_custom, inv_rate, pregenerated_tree, SEED, root_seq, "")
                        os.system(cur_cmd)
                        os.remove(prefix + ".full.treefile")
    log_file.close()

if __name__ == "__main__":
    main()
    # for _ in range(50):
    #     f = beautify_array(generate_f(True))
    #     print(f)

# iqtree2 --alisim generations/gtr_100_rnd_replica_0 -m gtr{0.01/0.1/0.21/0.26/0.05/0.18}+F{0.08/0.35/0.35/0.23}+G6+I{0.48} -t RANDOM{bd{0.1/0.05}/100} --seed 96734 --write-all --length 1140 -af fasta -rlen 0 0.001 0.01
# iqtree2 --alisim generations/gtr_1000_rnd_replica_0 -m gtr{0.11/0.32/0.26/0.94/0.23/0.25}+F{0.37/0.22/0.01/0.4}+G6+I{0.35} -t RANDOM{bd{0.1/0.05}/1000} --seed 96734 --write-all --length 1140 -af fasta -rlen 0 0.001 0.01
# iqtree2 --alisim generations/gtr_mam_rnd_replica_0 -m gtr{0.03/0.26/0.09/0.14/0.06/0.08}+F{0.26/0.38/0.16/0.21}+G6+I{0.04} -t /home/kpotoh/nemu-pipeline/data/alisim/external/mam_with_outgrp.nwk --seed 96734 --write-all --length 1140 -af fasta 
# iqtree2 --alisim generations/gtr_100_cytb_replica_0 -m gtr{0.3/0.16/0.19/0.26/0.01/0.04}+G6+I{0.35} -t RANDOM{bd{0.1/0.05}/100} --seed 96734 --write-all --root-seq /home/kpotoh/nemu-pipeline/data/alisim/external/human_cytb.fasta,human_cytb -af fasta -rlen 0 0.001 0.01
# iqtree2 --alisim generations/gtr_1000_cytb_replica_0 -m gtr{0.1/0.06/0.06/0.22/0.3/0.1}+G6+I{0.2} -t RANDOM{bd{0.1/0.05}/1000} --seed 96734 --write-all --root-seq /home/kpotoh/nemu-pipeline/data/alisim/external/human_cytb.fasta,human_cytb -af fasta -rlen 0 0.001 0.01
# iqtree2 --alisim generations/gtr_mam_cytb_replica_0 -m gtr{0.07/0.02/0.11/0.1/0.04/0.16}+G6+I{0.32} -t /home/kpotoh/nemu-pipeline/data/alisim/external/mam_with_outgrp.nwk --seed 96734 --write-all --root-seq /home/kpotoh/nemu-pipeline/data/alisim/external/human_cytb.fasta,human_cytb -af fasta 
# iqtree2 --alisim generations/12.12_100_rnd_replica_0 -m 12.12{-0.27/0.04/-0.72/-0.62/-0.08/-0.57/-0.82/-0.4/0.07/0.28/0.97}+F{0.28/0.4/0.02/0.3}+G6+I{0.09} -t RANDOM{bd{0.1/0.05}/100} --seed 96734 --write-all --length 1140 -af fasta -rlen 0 0.001 0.01
# iqtree2 --alisim generations/12.12_1000_rnd_replica_0 -m 12.12{-0.44/0.2/-0.13/0.37/0.8/-0.38/-0.25/-0.34/-0.52/0.45/0.56}+F{0.16/0.07/0.33/0.44}+G6+I{0.43} -t RANDOM{bd{0.1/0.05}/1000} --seed 96734 --write-all --length 1140 -af fasta -rlen 0 0.001 0.01
# iqtree2 --alisim generations/12.12_mam_rnd_replica_0 -m 12.12{0.28/0.65/-0.13/0.25/0.5/-0.17/0.38/0.73/0.01/-0.78/0.31}+F{0.15/0.21/0.24/0.39}+G6+I{0.22} -t /home/kpotoh/nemu-pipeline/data/alisim/external/mam_with_outgrp.nwk --seed 96734 --write-all --length 1140 -af fasta
# iqtree2 --alisim generations/12.12_100_cytb_replica_0 -m 12.12{-0.49/0.53/-0.13/-0.47/0.26/-0.41/-0.31/0.15/-0.27/0.63/-0.61}+G6+I{0.15} -t RANDOM{bd{0.1/0.05}/100} --seed 96734 --write-all --root-seq /home/kpotoh/nemu-pipeline/data/alisim/external/human_cytb.fasta,human_cytb -af fasta -rlen 0 0.001 0.01
# iqtree2 --alisim generations/12.12_1000_cytb_replica_0 -m 12.12{0.16/-0.97/0.06/0.14/-0.21/-0.03/0.93/-0.07/0.83/0.47/-0.42}+G6+I{0.3} -t RANDOM{bd{0.1/0.05}/1000} --seed 96734 --write-all --root-seq /home/kpotoh/nemu-pipeline/data/alisim/external/human_cytb.fasta,human_cytb -af fasta -rlen 0 0.001 0.01
# iqtree2 --alisim generations/12.12_mam_cytb_replica_0 -m 12.12{-0.16/0.15/0.92/0.39/-0.31/-0.45/-0.54/-0.43/-0.05/0.14/0.58}+G6+I{0.46} -t /home/kpotoh/nemu-pipeline/data/alisim/external/mam_with_outgrp.nwk --seed 96734 --write-all --root-seq /home/kpotoh/nemu-pipeline/data/alisim/external/human_cytb.fasta,human_cytb -af fasta