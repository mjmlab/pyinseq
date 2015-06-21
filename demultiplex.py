#!/usr/bin/env python
"""
Demultiplexes a FASTQ file into multiple files by 5' end barcode

Output file includes the input file name and the barcode:
e.g., file.fastq demultiplexes into file_CGAT.fastq, etc.


Generalize to read in the input file and the barcode file

Log the number of reads in the input file and the number in each barcode
(and the number and percent that do not belong to one of those barcodes)

"""

import gzip
from itertools import islice
import sys # for output logging
import os
import argparse

def demultiplex_fastq(fastq_file):
    if fastq_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open

    barcodes = { \
    'CGAT':'Input1', \
    'GCTA':'Input2', \
    'AGTC':'Output1', \
    'AAAA':'Output2'}

    with opener(fastq_file, 'r') as f:

        print('\n=== Demultiplexing FASTQ input file by 5\' barcode ===')

        (file_root, file_ext) = (os.path.splitext(fastq_file))

        # barcode length
        # in future generalize for other length barcodes
        # measure length of barcodes and make sure all have same length; else exit
        b_len = 4

        # fastq record
        # identifier =     record[0]
        # sequence =       record[1]
        # alt_identifier = record[2]
        # quality =        record[3]
        record = []
        for i,line in enumerate(f):
            record.append(line.rstrip('\n'))
            if i % 4 == 3:
                if record[1][0:b_len] in barcodes:
                    with open('{0}_{1}{2}'.format(file_root, record[1][0:b_len], file_ext), 'a') as fo:
                        fo.write('\n'.join(record) + '\n')
                record = []
            if (i+1) % 1E+5 == 0:
                if (i+1) % 1E+6 == 0:
                    print('\n=== Demultiplexing FASTQ input file by 5\' barcode ===')
                print('{0} records processed.'.format(i+1)) # index i starts at 0
        print('{0} total records processed.'.format(i+1))


def main():
    demultiplex_fastq('E689_400k_lines.fastq')

if __name__ == "__main__":
    main()
