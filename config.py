""" Configuration file for pyinseq.py package """

# Experiment name
# Used as a prefix for naming files and samples
# Alphanumeric characters only
experiment = ''

# Genome file in GenBank format (concatenated GenBank if multiple contigs)
gb = ''

# Barcodes (5' barcodes)
sample_file = ''

# Transposon end sequences
# In many cases the transposon sequence will be the same for both ends.
# This sequence is filtered from the reads prior to mapping.
# For mariner transposons do not include the terminal TA. Mariner hops into
# TA dinucleotides and adding the TA here will cause those nucleotides to
# be filtered prior to quality checking and mapping. Defaults are sequences for
# typical mariner ends.
transposon_left = 'ACAGGTTG'
transposon_right = 'ACAGGTTG'

# Path to Bowtie
bowtie_build = '/Users/markmandel/bowtie-1.1.1/bowtie-build'
bowtie = '/Users/markmandel/bowtie-1.1.1/bowtie'
