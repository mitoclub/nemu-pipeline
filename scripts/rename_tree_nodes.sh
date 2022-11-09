tree=$1  # path to newick tree

# ~/Downloads/RAxML_nodeLabelledRootedTree.nwk
# get internal node labels | drop ROOT | iteratively generate map-file

nw_labels -L $tree | grep -v ROOT | xargs -I{} echo -e "{}\tNode{}" > /tmp/map.tsv

if [ `grep -c NodeNode /tmp/map.tsv` -eq 0 ]; then
    nw_rename $tree /tmp/map.tsv
else
    cat $tree
fi
