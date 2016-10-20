# Documentation

## Running `pyinseq`

In general, place all input files in the `data/` directory. Each time the software is run an experiment directory is created in `results/` that contains the analysis files for that experiment.

`$ python run.py -i <reads> -s <sample_list> -g <genbank_reference> -e <experiment_name>`

`python` - Tested with Anaconda Python 3.5

`pyinseq.py` - main script

## Required arguments

`-i`  `--input` - Illumina reads file in FASTQ format

`-s`  `--samples` - sample list where each line is the `sample_name`tab`barcode(4bp)`

`-g`  `--genome` - concatenated genbank file for all contigs for the genome

`-e`  `--experiment` - all results files will be created in the subfolder of this name

## Optional arguments

`-d`  `--disruption` - fraction of gene disrupted (0.0 - 1.0)

`-d 1.0` is default (insertions anywhere in gene are counted)

`-d 0.9` includes only insertions in 5'-most 90% of each gene for scoring disruption.
