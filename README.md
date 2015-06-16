# pyinseq

**Developing Python software with the following goals**

1. Reproduce the functionality in [Andy Goodman's INSeq Perl pipeline](http://www.nature.com/nprot/journal/v6/n12/extref/nprot.2011.417-S2.zip) (15 MB .zip download) with additional improvements listed below.
2. Learn to be a better Python programmer.


**Implemented**

3. Read in FASTQ files instead of SCARF. (Note that David Cronin implemented this previously for us in Perl but in a complicated way.)
4. Demultiplex into separate files per barcode - will facilitate data deposition and analysis when samples from multiple experiments are run in the same Illumina run.

**Need to implement for functionality**

3. Run Bowtie (or Bowtie2?) seemlessly from within the software instead of requiring complicated index setup.
4. Finish the map to genes a la Goodman.
6. Summarize analysis results in a single table (for all chromosomes/plasmids; and for all samples under analysis).

**Wish list**

5. Include RNA genes.
5. Ability to analyze Left & Right transposon ends separately.
5. Count TA sites per gene and normalize genes per number of TA sites.
7. LOESS normalization to correct for more hits at origin and fewer hits at terminus.
8. Statistical analysis of results (method TBD).
9. Plotting of results (specific plots TBD); likely an analysis of essential and conditionally essential genes.
10. Detailed report per sample and for entire analysis (e.g., fraction of reads that mapped successfully; variation from other samples))

<img src="https://cloud.githubusercontent.com/assets/8669125/8175914/1dcad01c-13b8-11e5-8ceb-1b4f64a99f13.png" width="500">
