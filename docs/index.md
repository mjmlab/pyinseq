# Overview

This is a Python implementation of the methodology used by [Brooks et al. (2014)](http://www.ncbi.nlm.nih.gov/pubmed/25404340) to map transposon insertions in *Vibrio fischeri*. It follows closely from that described in [Goodman et al. (2011)](http://www.ncbi.nlm.nih.gov/pubmed/22094732), with exceptions as noted below.

The software is intended to be user-friendly, requiring only the FASTQ reads file, the GenBank file for the organism, and the list of samples and 5' DNA barcodes (barcodes can be omitted if the samples are already demultiplexed). Bowtie installation, transposon mapping, and gene mapping occur in an automated fashion.

Contents:

Core Documentation:
- [Installation](installation.md)
- [Documentation](documentation.md)

Additional Information:
- [Example data sets](example.md)   
- [Automated tests](tests.md)    
- [Roadmap](roadmap.md)   
- [Troubleshooting](troubleshooting.md)   
- [Comparison to other software](comparison.md)  

# License

BSD license as detailed in the [license file](../LICENSE.md).
