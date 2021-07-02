[![Build Status](https://travis-ci.com/mjmlab/pyinseq.svg?branch=master)](https://travis-ci.org/mjmlab/pyinseq)
![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)
![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)

<img src="img/pyinseq-logo.png" alt="drawing" width="200"/>

Lightweight python package to map transposon insertion sequencing (INSeq) data in
bacteria.

## IMPORTANT NOTE

We are merging in the `snakemake` & `conda` branch on July 2, 2021. `pyinseq` is not actually up on bioconda quite yet... It will be there very soon! In the meantime, this version is stable but cannot yet be installed via conda. You also may want to look at [v.0.2.1](https://github.com/mjmlab/pyinseq/releases/tag/v0.2.1). Contact us with any questions. We will remove this note once the package is available from bioconda.

## Quick start

This section is meant for users who know their way around terminal and `conda`. To use `pyinseq`, 
create a virtual environment with `python` 3.6 and install `pyinseq` using `conda`.

```bash
$ conda create -c bioconda -n pyinseq-py36 pyinseq
$ conda activate pyinseq-py36
```

Verify your installation with `--test`

```bash
(pyinseq-py36) $ pyinseq --test
```

Now you can run `pyinseq`!

```bash
$ pyinseq -i <input file> -s <sample file> -g <genbank file> -e <experiment name>
```

## Table of contents

* [Background](#Background)
* [Overview of command line operation](#Overview-of-command-line-operation)
* [Notes on output](#Notes-on-output)
* [User Guide](#User-guide)
  * [Input files description](#Input-files-description)
  * [General usage](#General-usage)
  * [Rerunning pyinseq](#Rerunning-pyinseq)
  * [Specialized tasks](#Specialized-tasks)
  * [Output description](#Output-files)
* [Installation](#Installation)
  * [Requirements](#Requirements)
  * [Using conda](#Using-conda-(recommended))
  * [Using pip](#Using-virtualenv-and-pip)
  * [Install from source code](#Install-from-source-code)
  * [Testing](#Testing)
* [FAQ](#FAQ)
* [License](#License)


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

- Generate genome files in **gff3** format.

`-d` / `--disruption`

- Five-prime fraction of gene (`0.0` - `1.0`) that must be disrupted for the hit to be counted in the summary_gene_table. Often insertions at the 3' end of a gene do not disrupt function so it may be of interest to run the pipeline with a disruption value of `0.8` or `0.9`.
  
`--min_count`

- Minimum number of reads per insertion site.

`--max_ratio`

- Maximum ratio of left:right or right:left reads per insertion site. 

`--barcode_length`

- Length of barcode index that is expected in samples.

`--transposon_seq`

- DNA sequence of transposon that flanks reads.

`-t` / `--threads`

- Number of cores to use for execution

`--additional_params`

- Additional parameters that will get passed to `snakemake`.


## Output description

| File | Description |
| --- | --- |
| `<experiment name>-config.yml` | Configuration file with run parameters |
| `results/summary_gene_table.txt` | summary for entire experiment |
| `results/<sample>_sites.txt` (for each sample) | Counts of each insertion in each sample |
| `results/<sample>_genes.txt` (for each sample) | Counts of each insertion mapped to genes |
| `results/<sample>_bowtie.txt` (for each sample) | Bowtie mapping results |
| `results/<sample>_trimmed.fastq` (for each sample) | Demultiplexed fastq reads trimmed for the chromosome sequence only |
| `results/log.txt` | text printed to console |
| `results/summary_log.txt` | summarized version of log output |
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

**summary_log.txt**
* Log file that will save printed output from `snakemake`. This includes the order of steps taken during the pyinseq execution.

**`<experiment name>`-config.yaml**
* Configuration file with run parameters.

## Background

`pyinseq` was inspired by the software published by [Goodman et al. (2011)](http://www.ncbi.nlm.nih.gov/pubmed/22094732). There are a number of differences, some of which are noted here. Most are fairly superficial in that they are intended to increase automation and reproducibility but do not materially affect the results. One exception is the first point, which does affect the resulting output.

1. As was conducted in [Brooks et al. (2014)](http://www.ncbi.nlm.nih.gov/pubmed/25404340) transposon site data are normalized across all replicons (chromosomes/plasmids/contigs) for calculation of a counts-per-million (cpm) per sample.

1. The user's file of reads is demultiplexed into separate .fastq files for each barcoded read.

1. The user provides a GenBank file, and `pyinseq` generates a fasta nucleotide file (.fna) and gene feature table (.ftt) that are formatted and named properly for `bowtie`.

1. The gene feature table (.ftt) is comparable to the protein feature table (.ptt) but it also includes RNA genes.

1. The `bowtie` software is installed through `conda`.

1. At the end of the analysis results are aggregated into a tab-delimited table and sample info is summarized.

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
additional_params: []
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
$ pyinseq -i reads.fastq -s samples.txt -g genome.gbk -e demo-run --threads 4 --additional_params --use-conda 
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
- [Python 3.6](https://www.python.org/downloads/release/python-3613/) (or 3.7)
- [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) (v1.3.0)

> `pyinseq` has not being tested on Windows operating systems, but as of Windows 10 
> there is support for terminals with *Ubuntu*

> Also note that `pyinseq` uses `bowtie` and **not bowtie2** which is a different software.

### Using conda (recommended)

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

#### Creating a virtual environment

A virtual environment is an isolated computational space where you can install dependencies and software without affecting the base operating system's configuration. 
We can use `conda` to create a virtual environment with python 3.6

```bash
$ conda create -n pyinseq python=3.6
```

> To use `python` 3.7, change `python=3.6` to `python=3.7`

To activate your environment:

```bash
$ conda activate pyinseq
```

You should see the name of your environment surrounded by parentheses in your terminal prompt.

```bash
(pyinseq) $ 
```
#### Installing `pyinseq` through `bioconda`

Now, using conda you can install `pyinseq` directly from the bioconda channel into your virtual environment.

```bash
(pyinseq) $ conda install -c bioconda pyinseq 
```

Verify that pyinseq installed correctly by running:

```bash
(pyinseq) $ pyinseq --help
2021-05-26 13:10 - INFO - pyinseq - Process command line arguments
usage: pyinseq [-h] [--get_default_config] [--config_format CONFIG_FORMAT]
               [-c CONFIG] [-t THREADS] [--additional_params ...] [-v]
               [-i INPUT] [-s SAMPLES] [-e EXPERIMENT] [-g GENOME]
               [-d DISRUPTION] [--min_count MIN_COUNT] [--max_ratio MAX_RATIO]
               [--barcode_length BARCODE_LENGTH]
               [--transposon_seq TRANSPOSON_SEQ] [--gff3]
               {demultiplex,genomeprep} ...
...
```
 Now you are ready to run `pyinseq`!

### Using `virtualenv` and `pip`

#### Install `virtualenv` and create an environment

If conda is not available, you can manually install [Python 3.6](https://www.python.org/downloads/release/python-3613/) (or 3.7) and use `pip` to install `virtualenv`. 

```bash
$ pip3 install virtualenv
```

Then use `virtualenv` to create a virtual environment. First, determine where your Python lives.

```bash
$ which python3
/Library/Frameworks/Python.framework/Versions/3.6/bin/python3
```

Then use this path to point to the Python interpreter that will be used in the virtual environment called `pyinseq`.

```bash
$ virtualenv -p /Library/Frameworks/Python.framework/Versions/3.6/bin/python3 ~/venvs/pyinseq
```

Activate this environment and install `pyinseq` using `pip`.

```bash
$ source ~/venvs/pyinseq/bin/activate
(pyinseq) $ pip install pyinseq
```

#### Install `bowtie` manually

Installation of `pyinseq` from `pip`/PyPi rather than conda will be missing [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) (v1.3.0). 
Download the binary release [here](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.0/) or use `curl`. For example, below I am using `curl` to download the zipfile into my computer.

```bash
# Download the zipfile with executables
(pyinseq) $ curl -L https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.0/bowtie-1.3.0-macos-x86_64.zip -O ~/bowtie-1.3.0-macos-x86_64.zip
(pyinseq) $ unzip bowtie-1.3.0-macos-x86_64.zip
```

Once you download the executables, you can add the `bin` folder to the `PATH` environment variable.

```bash
(pyinseq) $ export PATH=~/bowtie-1.3.0-macos-x86_64/:$PATH
```

Verify that `bowtie` executables are available on the terminal.

```bash
(pyinseq) $ bowtie --help
usage: bowtie [-h] [-b | -i] [--verbose] [--debug] [--large-index]
              [--index INDEX]
```

### Install from source code

If the above methods are not for you, you can directly clone the repository from github and install `pyinseq`

```bash
$ git clone https://github.com/mjmlab/pyinseq
$ cd pyinseq
$ pip install ./
```

Or just install directly from `github`.

```
pip install git+git://github.com/mjmlab/pyinseq
```

> Make sure that `bowtie` executables are available in your `PATH` variable. You can also follow this [section](#install-bowtie-manually) to install bowtie.

### Testing

You can test your installation of `pyinseq` by using the option `--test`.

```bash
$ pyinseq --test
```

If all tests pass then you are good to go!

## FAQ

### Why do I get errors in processing the GenBank file?

Ensure that the file is in GenBank and not RefSeq format.

### How do I uninstall pyinseq?

You can do this two ways: 
* uninstalling pyinseq from the virtual environment

```pip uninstall pyinseq``` or ```conda remove pyinseq```

* or completely remove the conda virtual environment.

``` conda env remove -n pyinseq```

### How can I notify of an issue with `pyinseq`

Please use the [GitHub Issues](https://github.com/mjmlab/pyinseq/issues).

## [License](LICENSE.md)

Pyinseq is an open-source software licensed under BSD-3.
