[![Build Status](https://travis-ci.org/mandel01/pyinseq.svg?branch=master)](https://travis-ci.org/mandel01/pyinseq)
![Python 3.5](https://img.shields.io/badge/python-3.5-blue.svg)

# pyinseq

Lightweight package to map transposon insertion sequencing (INSeq) data in
bacteria.

# Documentation

See the [pyinseq documentation](docs/index.md).

# Quickstart

## Installation

1. Download (or git clone) the `pyinseq` package using the link on this page.

1. Install the [Anaconda Python 3.5 download](https://www.continuum.io/downloads).

1. In the command line run the command below to install additional packages as
detailed in the [requirements.txt](requirements.txt) file.  
`pip install -r requirements.txt`

## Run the simple example dataset to check installation

`$ python run.py -i data/example/example01.fastq -s data/example/example01.txt -g data/example/ES114v2.gb -e example01`

### Required arguments

`-i`  `--input` - Illumina reads file in FASTQ format

`-s`  `--samples` - sample list where each line is the `sample_name`tab`barcode(4bp)`

`-g`  `--genome` - concatenated genbank file for all contigs for the genome

`-e`  `--experiment` - all results files will be created in the subfolder of this name

### Optional arguments

`-d`  `--disruption` - fraction of gene disrupted (0.0 - 1.0)

[Full documentation](docs/index.md)
