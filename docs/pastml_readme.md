## PastML

Used to reconstruct reference mutationals spectra

model HKY (maximal)

1. Reformat alignment for input. Genes separating

```bash
python mutspec/aln2pastml.py --aln data/example_nematoda/alignments_nematoda_clean --scheme data/example_nematoda/scheme_devilworm.nex --outdir data/example_nematoda/leaves
```

2. Repair scipy

There are error ([log](./pastml.log)) while minimizing some function during pastml run. Origin is scipy function in 
file `env_bio/lib/python3.9/site-packages/scipy/optimize/_numdiff.py`.

Lines 469-470 (`approx_derivative` func) replaced by

```python
    # if np.any((x0 < lb) | (x0 > ub)):
        # raise ValueError("`x0` violates bound constraints.")
    if np.any(((x0 < lb) | (x0 > ub)) & ~np.isclose(x0, lb) & ~np.isclose(x0, ub)):
        raise ValueError(
            "`x0` violates bound constraints. \nx0={}, \nlb={}, \nub={}, \n(x0 < lb)={}, \n(x0 > ub)={}, \nx0 type={}, \nlb type={}, \nub type={},"
            .format(x0, lb, ub, x0 < lb, x0 > ub, x0.dtype, lb.dtype, ub.dtype)
        )
```

3. Run pastml

```bash
pastml --prediction_method MPPA -m HKY -t iqtree_anc_tree.nwk -d mammals_nd1.tsv --work_dir pastml_states/ --html pastml_states/tree.html --threads 8

parallel  mkdir -p out_{/.} ';' pastml --prediction_method MPPA -m HKY -t {/.}.nwk -d {} --work_dir out_{/.} --html out_{/.}/tree.html --threads 8 ::: *.tsv
```

4. Reformat pastml output to usual states style

```bash
python mutspec/pastml2custom_format.py --model HKY --aln data/example_nematoda/alignments_nematoda_clean/ --outpath data/example_nematoda/genes_states.pastml_HKY.tsv data/example_nematoda/pastml_n_HKY/*
```