[![Build Status](https://travis-ci.com/mjmlab/pyinseq.svg?branch=master)](https://app.travis-ci.com/github/mjmlab/pyinseq)
![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)
![Python 3.8](https://img.shields.io/badge/python-3.8-blue.svg)

![logo](https://github.com/mjmlab/pyinseq/raw/master/img/pyinseq-logo-small.png)

Lightweight python package to map transposon insertion sequencing (INSeq) data in
bacteria.


## Quick start <!-- omit in TOC -->

This section is meant for users who know their way around terminal and `conda`. To use `pyinseq`, 
create a virtual environment with `python` 3.7 and install `pyinseq` using `conda`.

```bash
$ conda install -n base -c conda-forge mamba
$ conda create -n pyinseq-py37 python=3.7
$ conda activate pyinseq-py37
(pyinseq-py37) $ conda install -c bioconda bowtie
(pyinseq-py37) $ pip install pyinseq
```

Verify your installation with `--test`

```bash
(pyinseq-py37) $ pyinseq --test
```

Now you can run `pyinseq`!

```bash
(pyinseq-py37) $ pyinseq -i <input file> -s <sample file> -g <genbank file> -e <experiment name>
```

## Table of contents <!-- omit in TOC -->

- [Overview of Command Line Operation](#overview-of-command-line-operation)
- [Output description](#output-description)
- [Notes on output](#notes-on-output)
- [Background](#background)
- [User guide](#user-guide)
  - [Input files description](#input-files-description)
  - [General Usage](#general-usage)
  - [Specialized tasks](#specialized-tasks)
- [Installation](#installation)
  - [Requirements](#requirements)
  - [Installation in a conda environment (recommended)](#installation-in-a-conda-environment-recommended)
  - [Testing](#testing)
- [FAQ](#faq)
  - [Why do I get errors in processing the GenBank file?](#why-do-i-get-errors-in-processing-the-genbank-file)
  - [How do I uninstall pyinseq?](#how-do-i-uninstall-pyinseq)
  - [How can I notify of an issue with `pyinseq`](#how-can-i-notify-of-an-issue-with-pyinseq)
- [License](#license)


## Overview of Command Line Operation

Basic operation and a short description of the files are listed here. Below are detailed descriptions and links to example input files.

```bash
$ pyinseq -i <input file> -s <sample file> -g <genbank file> -e <experiment name>
```

`-i` / `--input`

- Illumina reads file in FASTQ or gzipped-FASTQ format.

`-s` / `--samples`

- Sample list where each line is the `sample_name` \<tab\> `barcode (4bp)`.


`-g` / `--genome`

- Concatenated GenBank File for all contigs for the genome.

`-e` / `--experiment`

- All results files will be created in a subfolder with this name.

#### Snakemake

`--get_default_config`

- Creates a default configuration file for running `pyinseq`

`--config_format`

- File format for configuration file (`yaml` or `json`)

`-c`/ `--config`

- Configuration file for running a `pyinseq` workflow. Every other argument will be ignored.

`--test`

- Runs pytest on installed `pyinseq` software.

#### Optional arguments

`--gff`

- Generate genome output file in **gff3** format.

`-d` / `--disruption`

- Five-prime fraction of gene (`0.0` - `1.0`) that must be disrupted for the hit to be counted in the summary_gene_table. Often insertions at the 3' end of a gene do not disrupt function so it may be of interest to run the pipeline with a disruption value of `0.8` or `0.9`. [default: 0.9]
  
`--min_count`

- Minimum number of reads per insertion site. [default: 3]

`--max_ratio`

- Maximum ratio of left:right or right:left reads per insertion site. [default: 10]

`--barcode_length`

- Length of barcode index that is expected in samples. [default: 4]

`--transposon_seq`

- DNA sequence of transposon that flanks reads. [default: 'ACAGGTTG']

`-t` / `--threads`

- Number of cores to use for execution [default: the CPU count of the computer]

`--snakemake_params`

- Additional parameters that will get passed to `snakemake`.


## Output description

| File | Description |
| --- | --- |
| `<experiment name>-config.yml` | Configuration file with run parameters |
| `results/summary_gene_table.txt` | summary for entire experiment (values in counts-per-million, cpm) |
| `results/<sample>_sites.txt` (for each sample) | Counts of each insertion in each sample |
| `results/<sample>_genes.txt` (for each sample) | Counts of each insertion mapped to genes |
| `results/<sample>_bowtie.txt` (for each sample) | Bowtie mapping results |
| `results/<sample>_trimmed.fastq` (for each sample) | Demultiplexed fastq reads trimmed for the chromosome sequence only |
| `results/log.txt` | text printed to console |
| `results/samples_info_yml` | basics stats for each sample | 
| `results/genome_lookup/genome.fna` | genome fasta nucleotide file |
| `results/genome_lookup/genome.ftt` | genome feature table | 
| bowtie indexes | index files created from genome by bowtie | 
| `results/raw_data/<sample>.fastq` (for each sample) | demultiplexed files for each sample/barcode |
| `results/raw_data/_other.fastq` | demultiplexed files for unrecognized barcodes |


## Notes on output

**T50**
* The minimum number of transposon insertion sites in the sample that account for at least 50% of the samples's reads. Used as a crude measure to detect bottlenecks, when comparing output T50 to the library input T50, or when comparing biological or technical replicates.

**samples_info.yml**
* Contains basic information of each sample in the experiment. 

**log.txt**
* Log file from pyinseq, including the messages from the `snakemake` execution. This includes the order of steps taken during the pyinseq execution.

**`<experiment name>`-config.yaml**
* Configuration file with run parameters.

## Background

`pyinseq` was inspired by the software published by [Goodman et al. (2011)](http://www.ncbi.nlm.nih.gov/pubmed/22094732). There are a number of differences, some of which are noted here. Most are fairly superficial in that they are intended to increase automation and reproducibility but do not materially affect the results. One exception is the first point, which does affect the resulting output.

1. As was conducted in [Brooks et al. (2014)](http://www.ncbi.nlm.nih.gov/pubmed/25404340) transposon site data are normalized across all replicons (chromosomes/plasmids/contigs) for calculation of a counts-per-million (cpm) per sample.

1. The user's file of reads is demultiplexed into separate .fastq files for each barcoded read.

1. The user provides a GenBank file, and `pyinseq` generates a fasta nucleotide file (.fna) and gene feature table (.ftt) that are formatted and named properly for `bowtie`.

1. The gene feature table (.ftt) is comparable to the protein feature table (.ptt) but it also includes RNA genes.

1. If `pyinseq` is installed from `conda`, then `bowtie` is installed from `conda` as a dependency.

1. At the end of the analysis, results are aggregated into a tab-delimited table and sample info is summarized.

1. `pyinseq` is written primarily in Python and uses `snakemake` as the workflow manager. You probably figured that out already.


## User guide

### Input files description

**Illumina sequencing reads [`-i`]** for the experiment. File can be uncompressed (.fastq or .fq) or gzip-compressed (.fastq.gz or .fq.gz).

[Example input file](https://github.com/mjmlab/pyinseq/blob/master/pyinseq/tests/data/input/example01.fastq)

```txt
@DGL9ZZQ1:720:C6YD0ACXX:2:1101:1246:2185 1:N:0:
GAAGCGACACCGACGACGGTAACAGGTTGGATGATAAGTCCCCGGTCTTCG
+
CCCFFFDDHHHHHCGHHHIJDHIIJJFHHJIJJJJIJJDHIJJHIAGIJJJ
...
```

**Sample file [`-s`]** describing the sample names and barcodes. Sample names should be restricted to letters, numbers, dash (-), and underscore (_), with a tab between the sample and the barcode in each row of a text file. It is recommended that the file be prepared in a text editor to ensure that additional hidden characters are not introduced. [BBEdit](http://www.barebones.com/products/bbedit/) is one option for Mac.

> Microsoft Excel can export tab-delimited files (.tsv), but do not use Microsoft Word for this purpose.

[Example sample file](https://github.com/mjmlab/pyinseq/blob/master/pyinseq/tests/data/input/example01.txt)

```txt
E001_01 GAAG
E001_02 CTTT
```

**GenBank file [`-g`]** listing the features and DNA sequence for the organism. If the organism has multiple chromosomes/contigs in the sequence the file they should be concatenated into a single file. Ensure that the double slash `//` at the end of the file remains to separate each contig.

> Files from NCBI GenBank often work where the corresponding files from NCBI RefSeq do not. Feel free to contact us with any questions here.

[Example GenBank file](https://github.com/mjmlab/pyinseq/blob/master/pyinseq/tests/data/input/ES114v2.gb)

```txt
LOCUS       CP000020             2897536 bp    DNA     circular BCT 02-APR-2008
DEFINITION  Vibrio fischeri ES114 chromosome I, complete sequence.
ACCESSION   CP000020
VERSION     CP000020.2  GI:171902228
KEYWORDS    .
SOURCE      Vibrio fischeri ES114
...
FEATURES             Location/Qualifiers
     source          1..2897536
                     /organism="Vibrio fischeri ES114"
                     /mol_type="genomic DNA"
                     /strain="ES114"
                     /db_xref="taxon:312309"
                     /chromosome="I"
     gene            complement(313..747)
                     /gene="mioC"
                     /locus_tag="VF_0001"
                     /old_locus_tag="VF0001"
     CDS             complement(313..747)
                     /dnas_title="FMN-binding protein MioC"
                     /gene="mioC"
                     /locus_tag="VF_0001"
                     /old_locus_tag="VF0001"
                     /codon_start=1
                     /transl_table=11
                     /product="FMN-binding protein MioC"
                     /protein_id="AAW84496.1"
                     /db_xref="GI:59478709"
                     /translation="MKKVSIITGSTLGGAEYVGDHLADLLEEMDFSTDIHNQPNLDDI
                     DIDSLWLLVVSTHGAGDYPDNIKPFIQQLESVTQPLSSVEFAVVAIGDSSYDTFCAAG
                     KSLQNTLKEHGAIEKYPLLEIDVTQNSIPEEPAELWLKQHIC"
...
ORIGIN
        1 aagatcactt aatatatata agatctttta aagagatctt ttattagatc tattatatag
       61 atcgtcgatc tctgtggata agtgataaat gatcaatagg atcatatact ttagatggat
      121 ccaaagttgt tatctttctt tgatcttcga tcggacagct tgaggacaaa agagttagtt
      181 atccacaagg ggggagggcg ttagatctta ttcaatggat aactataact tgatcactgg
      241 atctttctat agttatccac atagtaggta tcatctattt aataactttt atagatcgga
      301 caacactttt tattaacaaa tgtgttgttt tagccacaat tctgctggtt cttcagggat
      361 actattttga gttacatcta tctctaatag agggtacttt tcgatagcgc catgctcttt
      421 taaggtattt tgaagtgatt ttcctgctgc gcagaaagtg tcataacttg aatcaccgat
      481 agcaacaaca gcaaattcaa cgctcgataa tggctgagtc acgctttcta actgttgaat
      541 aaatggttta atgttgtcag ggtaatcacc agcaccgtga gttgatacaa cgagtaacca
      601 taagctatca atatcaatat catctaaatt tggttgatta tgaatatcgg tggaaaaatc
      661 catttcttct aataaatcag caagatggtc accaacatac tcagcaccac ctagagtgct
      721 tcctgtaata atagatactt ttttcatgaa tttatcctat aaaaatataa aaaatgggcc
      781 tacataggcc cattattaat cttattaata ttggttttat ttaccaatac agaatgaagt
...
//
LOCUS       CP000021             1330333 bp    DNA     circular BCT 02-APR-2008
DEFINITION  Vibrio fischeri ES114 chromosome II, complete sequence.
...
```

### General Usage

Pyinseq uses [snakemake](https://snakemake.readthedocs.io/en/stable/) to execute the following workflow: 

`pyinseq`

- Demultiplex a file of FASTQ reads.
- Write separate trimmed versions of the files (no barcodes, no transposon sequence).
- Map the trimmed reads to the genome.
- Quantify insertions per site and per ORF in the genome.
- Output summary gene table and report files describing the dataset.

The benefit of using `snakemake` is that it allows for parallel execution if more than one `thread` is provided.
The modularity of `snakemake` enables future feature implementation.  

#### Run Parameters (Config File)

`pyinseq` writes a configuration file (`<experiment name>-config.yaml`) into the working directory, which holds the parameters for the run. 
For the example command below, the `pyinseq-config.yaml` file will store the parameters used for this run.

```bash
$ pyinseq -i reads.fastq -s samples.txt -g genome.gbk -e pyinseq --threads 4
```

[Example Config File](https://github.com/mjmlab/pyinseq/blob/master/pyinseq/tests/data/input/pyinseq-config.yaml)

```bash
$ cat demo-run-config.yaml
snakemake_params: []
barcode_length: 4
command: pyinseq
config: false
config_format: yaml
disruption: 0.9
experiment: demo-run
genome: genome.gbk
get_default_config: false
gff3: false
input: reads.fastq.gz
max_ratio: 10
min_count: 3
samples: samples.txt
threads: 4
transposon_seq: ACAGGTTG
```

If you installed **conda** in your terminal, you can use the option `--use-conda` which will install bowtie during the executuion of the workflow. 

```bash
$ pyinseq -i reads.fastq -s samples.txt -g genome.gbk -e demo-run --threads 4 --snakemake_params --use-conda 
```

You can also get a default configuration file by using `--get_default_config` and modify it using a text-editor.

```bash
$ pyinseq --get_default_config
$ ls
default-config-pyinseq.yaml
```

> Make sure that all file paths in the configuration are correct


### Specialized tasks

These commands are useful when combining samples from multiple Illumina runs in a single `pyinseq` analysis. 
Both commands will also create a configuration file specific to each subcommand.

`pyinseq demultiplex`

- Demultiplex a file of Illumina reads.
- Writes separate trimmed versions of the files (no barcode, no transposon sequence) unless the optional `--notrim` flag is added.

```bash
$ pyinseq demultiplex -i <input file> -s <sample file> -e <experiment name>
```

`pyinseq genomeprep`

- Prepare the fasta nucleotide (.fna) and feature table (.ftt) files from a GenBank file.
- Is a good quick test to run on new GenBank files.
- Generates bowtie indexes unless the `--noindex` flag is added.
- Optionally create GFF3 file with the `-gff` flag (for use in other programs).

```bash
$ pyinseq genomeprep -g <genbank file> -e <experiment name>
```

> Similar to `pyinseq` main usage, you can run these commands using a configuration 
> file and pass additional `snakemake` options.



## Installation

### Requirements

`pyinseq` was written and tested using a MacOS (or Linux-based) operating system.

- MacOS or Linux-based operating systems.
- [Python 3.7](https://www.python.org/downloads/release/python-3713/)
- [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) (v1.0.0)
- [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

> `pyinseq` has not being tested on Windows operating systems, but as of Windows 10 
> there is support for terminals with *Ubuntu*

> Also note that `pyinseq` uses `bowtie`; **not bowtie2**, which is a different software.

### Installation in a conda environment (recommended)

Conda is a command-line package manager that can create virtual environments with 
all the necessary dependencies to run `pyinseq`. You can acquire conda by installing 
[Anaconda](https://docs.anaconda.com/anaconda/install/) and all the packages it 
brings, or by installing its lightweight version called 
[miniconda](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/macos.html). 
We recommend **miniconda** since it installs a minimal number packages, 
and is actually all you need to run `pyinseq`

Once conda is installed, you can verify it by running:

```bash
$ conda --help
```

#### Installing mamba

`mamba` is now required for `snakemake`. Install it system-wide with this command. The subsequent steps are shown using `conda` since that is more common for the
scientific audience, but can largely be accomplished in `mamba` if the user desires.

```
conda install -n base -c conda-forge mamba
```

#### Creating a virtual environment

A virtual environment is an isolated computational space where you can install dependencies and software without affecting the base operating system's configuration. 
We can use `conda` to create a virtual environment with python 3.7

```bash
$ conda create -n pyinseq-py37 python=3.7
```

To activate your environment:

```bash
$ conda activate pyinseq-py37
```

You should see the name of your environment surrounded by parentheses in your terminal prompt.

```bash
(pyinseq-py37) $ 
```
#### Install `bowtie`:

```bash
(pyinseq-py37) $ conda install -c bioconda bowtie
```
#### Install `pyinseq`:

Install the stable version:

```bash
(pyinseq-py37) $ pip install pyinseq
```

Or, install the most recent version from GitHub:

```
(pyinseq-py37) $ pip install git+https://github.com/mjmlab/pyinseq
```


Verify your installation with `--test`

```bash
(pyinseq-py37) $ pyinseq --test


Verify that pyinseq installed correctly by running:

```bash
(pyinseq-py37) $ pyinseq --help
2021-05-26 13:10 - INFO - pyinseq - Process command line arguments
usage: pyinseq [-h] [--get_default_config] [--config_format CONFIG_FORMAT]
               [-c CONFIG] [-t THREADS] [--snakemake_params ...] [-v]
               [-i INPUT] [-s SAMPLES] [-e EXPERIMENT] [-g GENOME]
               [-d DISRUPTION] [--min_count MIN_COUNT] [--max_ratio MAX_RATIO]
               [--barcode_length BARCODE_LENGTH]
               [--transposon_seq TRANSPOSON_SEQ] [--gff3]
               {demultiplex,genomeprep} ...
...
```
 Now you are ready to run `pyinseq`!

### Testing

You can test your installation of `pyinseq` by using the option `--test`.

```bash
(pyinseq-py37) $ pyinseq --test
```

If all tests pass then you are good to go!

## FAQ

### Why do I get errors in processing the GenBank file?

Ensure that the file is in GenBank and not RefSeq format.

### How do I uninstall pyinseq?

You can do this two ways: 
* uninstalling pyinseq from the virtual environment

```pip uninstall pyinseq```

* or completely remove the conda virtual environment.

``` conda env remove -n pyinseq-py37```

Note: will need to `conda deactivate` first to leave the environment.

### How can I notify of an issue with `pyinseq`

Please use the [GitHub Issues](https://github.com/mjmlab/pyinseq/issues).

## [License](LICENSE.md)

Pyinseq is an open-source software licensed under BSD-3.
