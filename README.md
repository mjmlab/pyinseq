# pyinseq

Lightweight package to map transposon insertion sequencing (INSeq) data in
bacteria.

This is currently a Python implementation of the methodology used by [Brooks et
al. (2014)](http://www.ncbi.nlm.nih.gov/pubmed/25404340) to map transposon
insertions in Vibrio fischeri. It follows closely from that described in
[Goodman et al. (2011)](http://www.ncbi.nlm.nih.gov/pubmed/22094732), except
that insertions are normalized across all replicons instead of within each
replicon.

The software is intended to be user-friendly, requiring only the FASTQ reads
file, the GenBank file for the organism, and the list of samples and 5' DNA
barcodes (barcodes can be omitted if the samples are already demultiplexed).
Bowtie installation, index formation, and gene mapping occurs in an automated
fashion.

# Installation

1. Download (or git clone) the `pyinseq` package using the link on this page.

1. Install the [Anaconda Python 3.5 download](https://www.continuum.io/downloads).

1. In the command line run the command below to install additional packages as
detailed in the [requirements.txt](requirements.txt) file.  
`pip install -r requirements.txt`

# Running the software

In general, place all input files in the `pyinseq/` directory. Each time the software is run an experiment directory is created that contains the analysis files for that experiment.

`$ python pyinseq.py -i <reads> -s <sample_list> -g <genbank_reference> -e <experiment_name>`

`python` - Tested with Anaconda Python 3.5

`pyinseq.py` - main script

Required arguments:

`-i`  `--input` - Illumina reads file in FASTQ format

`-s`  `--samples` - sample list where each line is the `sample_name`tab`barcode(4bp)`

`-g`  `--genome` - concatenated genbank file for all contigs for the genome

`-e`  `--experiment` - all results files will be created in the subfolder of this name

Optional arguments:

`-d`  `--disruption` - fraction of gene disrupted (0.0 - 1.0)

`-d 1.0` is default (insertions anywhere in gene are counted)

`-d 0.9` includes only insertions in 5'-most 90% of each gene for scoring disruption.

# Example data sets included

**example01 : Contrived data to test software with a limited dataset, including in this repo.**

`$ python pyinseq.py -i _exampleData/example01.fastq -s _exampleData/example01.txt -g _exampleData/ES114v2.gb -e example01`

**example02 : Large dataset of 100 million reads, dataset must be downloaded separately.**

[Example details](_exampleData/exampleData.md)

# Troubleshooting

Please use the GitHub Issues. This software is under active development.

### FAQ

1. Errors in the genome file: Use files from GenBank, not from RefSeq.

# License

BSD license as detailed in the [license file](LICENSE.md).
