#!/usr/bin/env python
"""
Demultiplexes a FASTQ file into multiple files by 5' end barcode

Output file includes the input file name and the barcode:
e.g., file.fastq demultiplexes into file_CGAT.fastq, etc.
Unrecognized barcodes get written to file_Other.fastq

Future:
Do barcode errors work to exit the program?
Flush the output files before appending to them.
Add argument parsing.
Output report should be in same order as input barcodes -- from Collections import OrderedDict ?
Include option to not write unassigned barcodes (would run faster).
"""

import gzip
import sys  # temporary - for command line arg
import os
import argparse

def barcodes_prep(samples):
    """
    Extract barcodes from an INSeq sample list and conduct basic quality checks

    Input samples is a tab-delimited file. On each line the sample identifier
    should be listed (no spaces), then after a tab the DNA barcode.
    Comment lines begin with '#' and are not read. E.g.:
    # Experiment E01 sample file by M. Mandel
    Input1<tab>CGAT
    Input2<tab>GCTA
    The function checks that all barcodes are the same length.
    The function returns a list of barcodes and sample names and the
    length of the barcodes.

    """

    # Extract barcode list from tab-delimited sample file.
    # Ensure uppercase and stripped
    print('\n===== Checking barcodes =====')
    barcode_list = {}
    for line in samples:
        if not line.startswith('#'):
            new_barcode = line.rstrip().split('\t')[1]
            if new_barcode in barcode_list:
                print('Error: redundant barcode {}'.format(new_barcode))
                exit(1)
            new_sample = line.rstrip().split('\t')[0]
            if new_sample in barcode_list.values():
                print('Error: redundant sample identifier {}'.format(new_sample))
                exit(1)
            barcode_list[new_barcode] = new_sample
            barcode_length = len(new_barcode)

    # Print barcodes and check for same length
    for b in sorted(barcode_list):
        print('{0}\t{1}'.format(b, barcode_list[b]))
        if len(b) != barcode_length:
            print('Error: barcodes are not the same length')
            exit(1)

    print('n={0} unique barcodes of length {1} nt'.format(len(barcode_list), barcode_length))

    return barcode_list.keys(), barcode_length



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

        file_prefix = 'Exp001'   # supply in command line arguments; check alphanumeric only
        if file_prefix:
            if not file_prefix.isalnum():
                print('Error: File prefix should be alphanumeric only. You entered {}'.format(file_prefix))
                exit(1)
            file_prefix += '_'
            print(file_prefix)

        # truncates the file with the 'w' option before writing data below
        # in cases where the output file already exists
        with open('{}assigned.fastq'.format(file_prefix), 'w') as fo:
            pass
        with open('{}assigned.fastq'.format(file_prefix), 'a') as fo:

            print('\n===== Assigning barcode and transposon information =====')

            # fastq record
            # identifier =     record[0]
            # sequence =       record[1]
            # alt_identifier = record[2]
            # quality =        record[3]
            record = []
            for i,line in enumerate(f):
                record.append(line.rstrip('\n'))
                if i % 4 == 3:

                    # barcode sequence in read
                    bc = record[1][0:b_len]
                    if bc in barcodes:
                        count_list[bc] += 1

                    # Tn location as number of nuceotides after the barcode
                    tn_loc = record[1].find('ACAGGTTG') - b_len # Tn location
                    tn_end = 'I'    # Left / Right / Identical
                    ta = 'N'    # Insertion at a TA dinucleotide
                    if (record[1][(tn_loc+b_len-2):(tn_loc+b_len)]) == "TA":
                        ta = 'Y'

                    print('{0}/barcode:{1}/Tn:{2}_{3}_{4}'.format(record[0],bc,tn_loc,ta,tn_end))
                    print(record[1][b_len:(tn_loc+b_len)])
                    print(record[2])
                    print(record[3][b_len:(tn_loc+b_len)])

                    fo.write('{0}/barcode:{1}/Tn:{2}_{3}_{4}\n{5}\n{6}\n{7}\n'.format(
                        record[0],bc,tn_loc,ta,tn_end,
                        record[1][b_len:(tn_loc+b_len)],
                        record[2],
                        record[3][b_len:(tn_loc+b_len)]))

                    record = []

    print(count_list)

            #print record



"""            if i % 4 == 3:
                # Write FASTQ record (4 lines) to file specific to its barcode
                if record[1][0:b_len] in barcodes:
                    with open('{0}_{1}{2}'.format(file_root, record[1][0:b_len], '.fastq'), 'a') as fo:
                        fo.write('\n'.join(record) + '\n')
                    # Count the read for the barcode
                    count_list[record[1][0:b_len]] += 1
                # Write unassigned barcodes to (root)_Other(extension) file
                else:
                    with open('{0}_{1}{2}'.format(file_root, 'Other', '.fastq'), 'a') as fo:
                        fo.write('\n'.join(record) + '\n')
                record = []
            if (i+1) % 4E+5 == 0:
                if (i+1) % 4E+6 == 0:
                    print('\n===== Demultiplexing FASTQ input file by 5\' barcode =====')
                # index i starts at 0
                # 4 lines per FASTQ record
                print('{0} records processed.'.format((i+1)/4))
        print('{0} total records processed.'.format((i+1)/4))

        print('\n===== Demultiplexing Summary =====')

        print('{0:,} out of {1:,} reads ({2:.1%}) were from a listed barcode'.format(sum(count_list.values()), (i+1)/4, float(sum(count_list.values()))/((i+1)/4)))
        print('\nbarcode\trecords\t%_from_list')
        for c in count_list:
            print('{0}\t{1:,}\t{2:.1%}'.format(c,count_list[c], float(count_list[c])/sum(count_list.values()), float(100*count_list[c])/((i+1)/4)))
        print('Other\t{0:,}\t'.format(((i+1)/4)-sum(count_list.values())))
"""

# ===== Start here ===== #

def main():
    fastq_file = sys.argv[1]
    sample_file = sys.argv[2]
    demultiplex_fastq(fastq_file, sample_file)

if __name__ == "__main__":
    main()
