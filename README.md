# scdeth manuscript figure generation and Snakemake pipeline

## Requirements
The [scdepth](https://github.com/gwlab-ca/scdepth) python package must be installed.

## Notes

### Running the Snakefile

The Snakefile contains the pipeline used to generate all analyses from the raw BAM, GTF, position, and probe inputs.

By default, samples without precomputed cache files are skipped unless the corresponding raw inputs have been unpacked locally. A subset of the raw inputs needed for cache generation can be downloaded from Zenodo. After downloading, extract any desired raw_*.tar archives into the root of this repository.

```
# Any subset of these archives may be used
tar -xvf raw_scrna.tar
tar -xvf raw_scrna_flex.tar
tar -xvf raw_visium.tar
tar -xvf raw_visium_hd.tar

snakemake --rerun-trigger mtime -j <jobs>
```

### Generating the figures

Running the Snakefile is not required to generate the manuscript figures. Only processed.tar is needed.

The Jupyter notebooks used to generate all figures and Table S1 are located in notebooks/. Run notebooks/TableS1.ipynb first; the remaining notebooks can then be run in any order.

The `processed.tar` archive, available from Zenodo, contains all data required to run the notebooks. Extract it into the root of this repository:

```
tar -xvf processed.tar
jupyter lab
# Run notebooks/TableS1.ipynb first, then the remaining notebooks in any order
```
