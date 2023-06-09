#!/bin/bash

tree=$1
my_tmp_map=`mktemp`
my_tmp_tree=`mktemp`

nw_labels -L $tree | xargs -I{} echo -e "{}\tNode{}" > $my_tmp_map
nw_rename $tree $my_tmp_map > $my_tmp_tree
cat $my_tmp_tree > $tree
rm $my_tmp_tree $my_tmp_map
