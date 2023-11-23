## Containers

1. nemu-env.def - last and used version of container for nemu pipeline execution
2. bioinfo.def - different tools for phylogenetic analyses including nemu pipeline execution

## Build

```bash
singularity build image_nemu.sif nemu-env.def
```

## Run

```bash
singularity shell image_pipeline_latest.sif
# OR
./image_pipeline_latest.sif
```
