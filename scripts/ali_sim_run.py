import os
import sys
import random

N = 10
THREADS = 4
SEED = random.randint(0, 100000)
MAMMALS_TREE = "/home/kpotoh/dolphin/data/alisim/external/mam.nwk"
ROOT_SEQ = "/home/kpotoh/dolphin/data/alisim/external/human_cytb.fasta"

cmd = "iqtree2 --alisim {} -m {}{}+F{}+G6+I{} -t {} --seed {} --write-all -nt {} {} -af fasta"


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
        raise NotImplementedError
    
    return r


# GTR
mdl = "gtr"

for label in ["rnd", "cytb"]:
    root_seq = "" if label == "rnd" else f"--root-seq {ROOT_SEQ},human_cytb"

    # for raw_tree in ["mam", 100, 1000]:
    for raw_tree in ["mam"]:
        if isinstance(raw_tree, int):
            tree = "RANDOM{bd{0.1/0.05}/" + str(raw_tree) + "}"
        elif isinstance(raw_tree, str):
            tree = MAMMALS_TREE
        else:
            raise NotImplementedError

        for i in range(N):
            # mdl_args = "{0.5/0.6/0.9/0.2/0.1/0.4/0.7/0.8/0.3/0.15/0.65}"
            mdl_args = beautify_array(generate_rates(6, True))
            f_custom = beautify_array(generate_f())
            inv_rate = beautify_array([random.random() * 0.5])  # A proportion of invariable sites
            cur_cmd = cmd.format(f"generations/{mdl}_{raw_tree}_{label}_replica_{i}", mdl, mdl_args, f_custom, inv_rate, tree, SEED, THREADS, root_seq)
            print(cur_cmd)
            os.system(cur_cmd)