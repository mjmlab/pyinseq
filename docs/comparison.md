# Comparison with other software

`pyinseq` was inspired by the software published by [Goodman et al. (2011)](http://www.ncbi.nlm.nih.gov/pubmed/22094732). There are a number of differences, many of which are noted here. Most are fairly superficial in that they are intended to increase automation and reproducibility but do not materially affect the results. One exception is the first point, which does affect the resulting output.

1. As was conducted in [Brooks et al. (2014)](http://www.ncbi.nlm.nih.gov/pubmed/25404340) transposon site data are normalized across all replicons (chromosomes/plasmids/contigs) for calculation of a counts-per-million (cpm) per sample.

1. The user provides a GenBank file, and `pyinseq` generates a fasta nucleotide file (.fna) and gene feature table (.ftt) that are formatted and named properly for `bowtie`.

1. The user's file of reads is demultiplexed into separate .fastq.gz files for each barcoded read.

1. The gene feature table (.ftt) is comparable to the protein feature table (.ptt) but it also includes RNA genes.

1. The `bowtie` software is packaged with `pyinseq` (a link to its own license is provided on the [LICENSE.md page](../LICENSE.md).

1. At the end of the analysis results are aggregated into a tab-delimited table.

1. `pyinseq` is written primarily in Python. You probably figured that out already.
