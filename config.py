""" Configuration file for pyinseq.py package """

# TRANSPOSON END SEQUENCES
# Typical mariner transposon ends for INSeq are the defaults (same left/right ends)
# This sequence is filtered from the reads prior to mapping
# Mariner transposons inserts at TA dinucleotides; do not include the corresponding
# TA from the tranpsoson here
transposonLeft = 'ACAGGTTG'
transposonRight = 'ACAGGTTG'

# PATH TO BOWTIE
bowtie = '/Users/markmandel/bowtie-1.1.1/bowtie'
bowtieBuild = '/Users/markmandel/bowtie-1.1.1/bowtie-build'
