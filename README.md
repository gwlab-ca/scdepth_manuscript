# scdeth manuscript figure generation and Snakemake pipeline

## Requirements
The [scdepth](https://github.com/gwlab-ca/scdepth) python package must be installed.

## Notes

 * The Snakefile contains the pipeline that was used to generate all of the analyses from the raw bam/gtf/position/probe data
 * The Snakefile is setup to ignore generating cached tags.  Raw cached tag data can be downloaded from [Zenodo](https://dx.doi.org/10.5281/zenodo.15518941). The `raw_*.tar` files must be extracted into the root of this repository
 * The jupyter notebooks to generate all of the Figures and Table S1 are in [notebooks](./notebooks)
 * The Table S1.ipynb notebook should be run first
 * The `processed.tar` file can also be downloaded from [Zenodo](https://dx.doi.org/10.5281/zenodo.15518941). It contains all the data necessary to run the notebooks
 * The `processed.tar` file must be extracted in the root of this repository
