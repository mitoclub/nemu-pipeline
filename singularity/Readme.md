## containers

1. nemu.def - last and used version of container while isolated pipeline execution
2. p.def - python dependencies
3. pipeline-2.5.def - legacy

### def files for fasta rebuilding

1. appendix.def - add pyvolve and mutspec-utils (legacy)
2. appendix2.def - add pymutspec, macse 2.06, mafft


## build
```bash
singularity build image_pipeline_latest.sif nemu.def
singularity build image_pipeline_latest.sif appendix2.def
```