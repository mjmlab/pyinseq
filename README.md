# pyinseq

## Description

Python package to map transposon insertion sequencing (INSeq) data in bacteria

# Current status

I've been using this in-house but there are still some bugs. Many steps are commented out.

# Project purpose

Reproduce the functionality in [Andy Goodman's INSeq Perl pipeline](http://www.nature.com/nprot/journal/v6/n12/extref/nprot.2011.417-S2.zip) (15 MB .zip download) in Python to serve as a platform for additional functionality for normalization, plotting, sample management, and data analysis.

Start at `pyinseq.py` and that calls the other scripts/modules as needed. Lines commented out at the far left need to be uncommented, then checked with a full run with test data.

# Implemented

1. Read in FASTQ files instead of SCARF. (Note that David Cronin implemented this previously for us in Perl but in a manner that is now yielding multiple errors.)
2. Demultiplex into separate files per barcode - will facilitate data deposition and analysis when samples from multiple experiments are run in the same Illumina run.
  - *This has since been removed from the main branch for speed but can be added back in as needed (perhaps with optimized software). It was not necessary to analyze samples from only a single run but will be necessary to analyze samples from separate runs.*
3. Run bowtie seemlessly from within the software instead of requiring complicated index setup. Improves upon the previous software by reading a GenBank file, generating appropriately names files for Bowtie, calling bowtie-build to build indexes, and then calling bowtie for mapping. Additionally the .ftt files generated mimic the .ptt file format but additionally include RNA genes.
4. Normalization. Do CPM normalization to mimic previous analysis. [**check math here**]
4. Map number of hits per gene.
5. Summarize analysis results in a single table (for all chromosomes/plasmids; and for all samples under analysis).

See the [Roadmap file](roadmap.md) for info on specific items to be implemented, including bugs in this version.

# Next steps/wish list

5. Add back in demultiplexing and general method to organize samples in folders.
6. Ability to analyze Left & Right transposon ends separately.
  - Some code for this is already in to note L/R side on the identifier line but it has not been used.
7. Count TA sites per gene and normalize gene size per number of TA sites.
7. LOESS normalization to correct for more hits at origin and fewer hits at terminus.
8. Statistical analysis of results (method TBD; consider DESeq2 and other RNA-seq methods).
9. Plotting of results (specific plots TBD); likely an analysis of essential and conditionally essential genes. Consider assigning samples.
10. Detailed reporting and logging (e.g., fraction of reads that mapped successfully; variation from other samples))
11. Overall optimization

# License

Want people to be able to adapt it to their needs with attribution here as appropriate. So I should find an appropriate license that captures that...
