#!/usr/bin/env python
"""
Demultiplexes a FASTQ file into multiple files by 5' end barcode

Output file includes the input file name and the barcode:
e.g., file.fastq demultiplexes into file_CGAT.fastq, etc.


Future:
Add argument parsing.
Output report should be in same order as input barcodes.

NOTE ERROR IN CALCULATIONS -- LISTS 40 RECORDS WHEN SHOULD LIST 10 (4x) DUE TO FASTQ LINE COUNTING.

"""

import gzip
import os
import argparse

def barcodes_prep(samples):
    """
    Extract barcodes from an INSeq sample list and conduct basic quality checks

    Input samples is a tab-delimited file. On each line the 5' barcode should
    be listed, then after a tab the description of the sample. E.g.:
    CGAT   Input1
    GCTA   Input2
    The function checks that all barcodes are the same length.
    The function returns a list of barcodes (no sample names) and the
    length of the barcodes.

    """

    # Extract barcode list from tab-delimited sample file.
    # Ensure uppercase and stripped
    barcode_list = []
    for line in samples:
        barcode_raw = (line.rstrip().split('\t'))[0]
        barcode_list.append(barcode_raw.upper().rstrip())

    barcode_length = len(barcode_list[0])

    print('\n===== Checking barcodes =====')

    # Print barcodes and check for same length
    for b in barcode_list:
        print(b)
        if len(b) != barcode_length:
            print('Error: barcodes are not the same length')
            exit() # How to really do error handling???

    if len(barcode_list) != len(set(barcode_list)):
        print('Error: non-unique barcodes in samples list')
        exit() # How to really do error handling???

    print('n={0} unique barcodes of same length ({1} nt)'.format(len(barcode_list), barcode_length))

    return barcode_list, barcode_length

def demultiplex_fastq(fastq_file, sample_file):
    if fastq_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open

    # barcodes = list of barcodes (sequences only)
    # b_len = barcode length
    with open(sample_file, 'r') as input_sample_file:
        barcodes, b_len = barcodes_prep(input_sample_file)

    # initialize to count reads per barcode
    count_list = {}
    for b in barcodes:
        count_list[b] = 0

    with opener(fastq_file, 'r') as f:

        print('\n===== Demultiplexing FASTQ input file by 5\' barcode =====')

        (file_root, file_ext) = (os.path.splitext(fastq_file))

        # fastq record
        # identifier =     record[0]
        # sequence =       record[1]
        # alt_identifier = record[2]
        # quality =        record[3]
        record = []
        for i,line in enumerate(f):
            record.append(line.rstrip('\n'))
            if i % 4 == 3:
                # Write FASTQ record (4 lines) to file specific to its barcode
                if record[1][0:b_len] in barcodes:
                    with open('{0}_{1}{2}'.format(file_root, record[1][0:b_len], file_ext), 'a') as fo:
                        fo.write('\n'.join(record) + '\n')
                    # Count the read for the barcode
                    count_list[record[1][0:b_len]] += 1
                record = []
            if (i+1) % 1E+5 == 0:
                if (i+1) % 1E+6 == 0:
                    print('\n===== Demultiplexing FASTQ input file by 5\' barcode =====')
                print('{0} records processed.'.format(i+1)) # index i starts at 0
        print('{0} total records processed.'.format(i+1))

        print('\n===== Demultiplexing Summary =====')

        print('{0:,} out of {1:,} reads ({2:.1%}) were from a listed barcode'.format(sum(count_list.values()), i+1, float(sum(count_list.values()))/(i+1)))
        print('\nbarcode\trecords\t%_from_list')
        for c in count_list:
            print('{0}\t{1:,}\t{2:.1%}'.format(c,count_list[c], float(count_list[c])/sum(count_list.values()), float(100*count_list[c])/(i+1)))
        print('Other\t{0:,}\t'.format((i+1)-sum(count_list.values())))

# ===== Start here ===== #

def main():
    fastq_file = 'E689_400k_lines.fastq'
    sample_file = 'samples.txt'
    demultiplex_fastq(fastq_file, sample_file)

if __name__ == "__main__":
    main()
