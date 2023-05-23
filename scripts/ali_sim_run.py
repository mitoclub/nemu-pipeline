import os
import sys
import random

N = 10
seed = random.randint(0, 100000)

cmd = "iqtree2 --alisim {} -m {}{}+F{}+I{} -t {} --seed {} --write-all"

# mdl = "10.12"
mdl = "gtr"
tree = "tree.nwk"


def beautify_array(lst):
    lst = [str(round(x, 2)) for x in lst]
    return "{" + "/".join(lst) + "}"


def generate_f() -> list:
    a = [random.random() for _ in range(4)]
    s = sum(a)
    a = [x / s for x in a]
    return a

def generate_rates(n=12, rev=True):
    a = [random. for _ in range(n)]


for i in range(N):
    mdl_args = "{0.5/0.6/0.9/0.2/0.1/0.4/0.7/0.8/0.3/0.15/0.65}"
    f_custom = beautify_array(generate_f())
    inv_rate = beautify_array([random.random() * 0.4])  # A proportion of invariable sites
    cur_cmd = cmd.format(f"aln_{i}", mdl, mdl_args, f_custom, inv_rate, tree, seed)
    print(cur_cmd)