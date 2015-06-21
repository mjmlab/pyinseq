""" Configuration file for pyinseq.py package

Edit each of the variables below as appropriate for current sample.
XXXXX IS THERE A MORE USER-FRIENDLY WAY TO POPULATION THESE FIELDS? XXXXX """



""" Barcodes for this sample

5' 4-nt barcodes """

barcodes = { \
'GTAC':'Input1','AAAA':'Input2', \
'TATA':'Input3','GTCA':'Input4', \
'TTAA':'Output5','AACC':'Output6', \
'GCTA':'Output7','AGTC':'Output8'}


""" Mariner transposon end sequences

In many cases the transposon sequence will be the same for both ends
XXXXX Insert comment about TA site """

transposon_left = ''
transposon_right = ''



""" Genome files """

fasta = ''
ptt = ''
rpt = ''
