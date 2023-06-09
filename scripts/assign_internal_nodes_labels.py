# USAGE: script.py TREE [TREE_LABELED]

import sys
from Bio import Phylo

assert len(sys.argv) > 1

intree  = sys.argv[1]
outtree = sys.argv[2] if len(sys.argv) == 3 else sys.stdout

# read the Newick file
tree = Phylo.read(intree, 'newick')

# traverse the tree and add internal node labels
for i, node in enumerate(tree.get_nonterminals(), 1):
    node.name = 'Node' + str(i)

# write the updated tree to a new file
Phylo.write(tree, outtree, 'newick')

