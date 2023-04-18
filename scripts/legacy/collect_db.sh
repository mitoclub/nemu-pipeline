

species=("homo" "gorilla" "mus" "lemur" "pan")
gene="ATP6"

script=/export/src/dolphin/scripts/geneUploader.pl
dir=/export/data/base

for i in ${!species[@]}; do
    name=${species[i]}
    singularity exec --bind /mnt/data/export:/export src/image_pipeline-2.5.sif perl $script $name $gene $dir
done
