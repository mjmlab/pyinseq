# pyinseq

Developing a Python software package with the following goals:

1. Learn to be a better Python programmer.
2. Reproduce the functionality in [Andy Goodman's INSeq Perl pipeline](http://www.nature.com/nprot/journal/v6/n12/extref/nprot.2011.417-S2.zip) (15 MB .zip download) with additional improvements listed below.
3. Run Bowtie (or Bowtie2?) seemlessly from within the software instead of requiring complicated index setup.
4. Demultiplex into separate files per barcode - will facilitate data deposition and analysis when samples from multiple experiments are run in the same Illumina run.
5. Summarize analysis results in a single table (for all chromosomes/plasmids; and for all samples under analysis).
6. LOESS normalization to correct for more hits at origin and fewer hits at terminus.
7. Statistical analysis of results (method TBD).
8. Plotting of results (specific plots TBD); likely an analysis of essential and conditionally essential genes.
9. Detailed report per sample and for entire analysis (e.g., fraction of reads that mapped successfully; variation from other samples)


