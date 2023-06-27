import os
import sys
import random

from ete3 import PhyloTree

MAMMALS_TREE = "/home/kpotoh/dolphin/data/alisim/external/mam_with_outgrp.nwk"
ROOT_SEQ = "/home/kpotoh/dolphin/data/alisim/external/human_cytb.fasta"
ROOT_SEQ_NAME = "human_cytb"
LENGTH = 1140
SEED = random.randint(0, 100000)

debug_mode = True
OUTDIR = "generations_mam"
N = 5


cmd = "iqtree2 --alisim {} -m {}{}{}+G6+I{} -t {} --seed {} --write-all {} -af fasta {} -redo"


if not os.path.exists(OUTDIR):
    os.mkdir(OUTDIR)
    print(f"{OUTDIR} created", file=sys.stderr)


def beautify_array(lst):
    lst = [str(round(x, 2)) for x in lst]
    return "{" + "/".join(lst) + "}"


def generate_f() -> list:
    a = [random.random() for _ in range(4)]
    s = sum(a)
    a = [x / s for x in a]
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


# mdl = "gtr"
for mdl in ["gtr", "12.12"]:
    for gene in ["rnd", "cytb"]:
        root_seq = f"--length {LENGTH}" if gene == "rnd" else f"--root-seq {ROOT_SEQ},{ROOT_SEQ_NAME}"

        for raw_tree in [100, 1000]:
        # for raw_tree in ["mam"]:
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
                f_custom = "+F" + beautify_array(generate_f()) if gene == "rnd" else ""
                inv_rate = beautify_array([random.random() * 0.5])  # A proportion of invariable sites
                cur_cmd = cmd.format(prefix, mdl, mdl_params, f_custom, inv_rate, tree, SEED, root_seq, rlen_args)
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
