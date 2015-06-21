""" Configuration file for pyinseq.py package """

# Barcodes (5' barcodes)
barcodes = { \
'CGAT':'Input1', \
'GCTA':'Input2', \
'AGTC':'Output1', \
'AAAA':'Output2'}

# Mariner transposon end sequences
# In many cases the transposon sequence will be the same for both ends
# XXXXX Insert comment about TA site
transposon_left = ''
transposon_right = ''

# Genome files. Files in same order for all sections below (e.g., Ch1, plasmid)
fasta = ''
ptt = ''
rnt = ''

# Path to Bowtie
bowtie = '/Users/markmandel/bowtie-1.1.1/bowtie'
bowtie-build = '/Users/markmandel/bowtie-1.1.1/bowtie-build'
