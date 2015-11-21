# pyinseq

# Description

Python package to map transposon insertion sequencing (INSeq) data in bacteria

**Warning to biologists...this is close, but not quite ready for prime time!**

# Development Priorities

Move toward a modular approach. General picture to split into three steps:

![modular approach](https://cloud.githubusercontent.com/assets/8669125/10409855/18925d7a-6ef5-11e5-9304-9f24eb868b80.png)

Update 11/21/15 - Modular approach in place.  
Step 1 could be faster but works.  
Step 2 works. Need to add (1) package bowtie together so separate installation not needed, (2) LOESS normalization, (3) filter by 5'-3' position, (4) filter by number of reads.  
Step 3. No code here yet.  
Additional reporting/logging needed throughout.  

# Installation

**Dependencies**

Recommend the [Anaconda Python 3.5 download](https://www.continuum.io/downloads). Step 3 steps will use `pandas` and other functions that Anaconda makes easy to install.  

[bowtie](http://bowtie-bio.sourceforge.net/index.shtml) - Latest version tested here is bowtie 1.1.1.

`screed` - `pip install screed`

**Download**

Download and open the zip file for the pyinseq package (or git clone).

**Config**

Edit the `config.py` file with the path to your bowtie installation and any other relevant changes.

# Running the software

In general, place all input files in the `pyinseq/` directory. Each time the software is run an experiment directory is created that contains the analysis files for that experiment.

`$ python pyinseq.py -i reads.fastq -s samples.txt -g genome.gb -e exp01`

`python` - Tested with Anaconda Python 3.5

`pyinseq.py` - main script

Required arguments:

`-i`  `--input` - concatenated genbank file for all contigs for the genome

`-e`  `--experiment` - all results files will be created in the subfolder of this name

`-i`  `--input` - Illumina reads file in FASTQ format

`-s`  `--samples` - sample list where each line is the `sample_name`tab`barcode(4bp)`

Optional arguments:

`-d`  `--disruption` - fraction of gene disrupted (0.0 - 1.0). **NOT CURRENTLY IN PLACE**

`-d 1.0` is default (insertions anywhere in gene are counted)

`-d 0.9` includes only insertions in 5'-most 90% of each gene for scoring disruption.

# Example data sets included

**example01 : Contrived data to test software with a limited dataset**

`$ python pyinseq.py -i _exampleData/example01.fastq -s _exampleData/example01.txt -g _exampleData/ES114v2.gb -e example01`

**example02 : large dataset of 100 million reads**

[Example details and additional examples](_exampleData/exampleData.md)

# Project purpose

Reproduce the functionality in [Andy Goodman's INSeq Perl pipeline](http://www.nature.com/nprot/journal/v6/n12/extref/nprot.2011.417-S2.zip) (15 MB .zip download) in Python to serve as a platform for additional functionality for normalization, plotting, sample management, and data analysis.

Immediate next steps will include multiple quality checks on the run and distinguishing the directionality of the transposon for a new transposon being built.


# Implemented

1. Read in FASTQ files instead of SCARF. (Note that David Cronin implemented this previously for us in Perl but in a manner that is now yielding multiple errors.)
2. Demultiplex into separate files per barcode - will facilitate data deposition and analysis when samples from multiple experiments are run in the same Illumina run.
3. Run bowtie seemlessly from within the software instead of requiring complicated index setup. Improves upon the previous software by reading a GenBank file, generating appropriately names files for Bowtie, calling bowtie-build to build indexes, and then calling bowtie for mapping. Additionally the .ftt files generated mimic the .ptt file format but additionally include RNA genes.
4. Normalization. Do CPM normalization to mimic previous analysis.
4. Map number of hits per gene.
5. Summarize analysis results in a single table (for all chromosomes/plasmids; and for all samples under analysis).

# Next steps/wish list

6. Ability to analyze Left & Right transposon ends separately.
  - Some code for this is already in to note L/R side on the identifier line but it has not been used.
7. Count TA sites per gene and normalize gene size per number of TA sites.
7. LOESS normalization to correct for more hits at origin and fewer hits at terminus.
8. Statistical analysis of results (method TBD; consider DESeq2 and other RNA-seq methods).
9. Plotting of results (specific plots TBD); likely an analysis of essential and conditionally essential genes. Consider assigning samples.
10. Detailed reporting and logging (e.g., fraction of reads that mapped successfully; variation from other samples))
11. Overall optimization

Current high-priority issues/bugs in Github issues. Slightly out-of-date thoughts and ideas are listed in the [Roadmap file](roadmap.md).

# License

BSD license as detailed in the [license file](LICENSE.md).
