#!/usr/bin/env python
"""
Trim a FASTQ file to removes barcode and transposon sequence.

Removes barcode (need to implement).
Removes tranpsoson sequence (need to implement).

Output file adds the modifier _trim:
e.g., file_CGAT.fastq becomes file_CGAT_trim.fastq

Future:
Add argument parsing.
Flush the output files before appending to them.
Options for different left/right transposon sequences.
Add handling for no transposon sequence present (now will fail at Bowtie)
Count number in which transposon sequence was found
Filter for transposon sequence TA site required?

"""

import gzip
import os
import argparse

def trim_fastq(fastq_file, barcode, tn_sequence):
    if fastq_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open

    #initialize to count number with transposon sequence
    count_list = {}

    with opener(fastq_file, 'r') as f:

        print('\n===== Trimming 5\' barcode and 3\' tn sequence =====\nTrimming file: {0}'.format(fastq_file))

        (file_root, file_ext) = (os.path.splitext(fastq_file))

        # fastq record
        # identifier =     record[0]
        # sequence =       record[1]
        # alt_identifier = record[2]
        # quality =        record[3]
        record = []

        for i,line in enumerate(f):
            if i % 2 == 0: # identifier or alt_identifier
                record.append(line.rstrip('\n'))
            if i % 4 == 1: # sequence
                # trim barcode and leading 'N' if present
                read_trimmed_bc = line.rstrip('\n')[len(barcode):].lstrip('N')
                # trim transposon sequence
                tn_site = read_trimmed_bc.find(tn_sequence)
                read_trimmed = read_trimmed_bc[:tn_site].rstrip('\n')
                record.append(read_trimmed)
            if i % 4 == 3: # quality
                # trim the quality scores
                # keep the values corresponding to the remaining bases
                start = (len(line.rstrip('\n')) - len(read_trimmed_bc))
                record.append(line[start : (start + tn_site)])

            # Write FASTQ record (4 lines) to file
            if i % 4 == 3:
                with open('{0}_{1}{2}'.format(file_root, 'trim', file_ext), 'a') as fo:
                    fo.write('\n'.join(record) + '\n')
                record = []
            if (i+1) % 4E+5 == 0:
                if (i+1) % 4E+6 == 0:
                    print('\n===== Trimming 5\' barcode and 3\' tn sequence =====\nTrimming file: {0}'.format(fastq_file))
                print('{0:,} records processed.'.format((i+1)/4)) # index i starts at 0
        print('{0:,} total records processed.'.format((i+1)/4))

"""        print('\n===== Demultiplexing Summary =====')

        print('{0:,} out of {1:,} reads ({2:.1%}) were from a listed barcode'.format(sum(count_list.values()), (i+1)/4, float(sum(count_list.values()))/((i+1)/4)))
        print('\nbarcode\trecords\t%_from_list')
        for c in count_list:
            print('{0}\t{1:,}\t{2:.1%}'.format(c,count_list[c], float(count_list[c])/sum(count_list.values()), float(100*count_list[c])/((i+1)/4)))
        print('Other\t{0:,}\t'.format(((i+1)/4)-sum(count_list.values()))) """

# ===== Start here ===== #

def main():
    fastq_file = 'E689_400k_lines_AAAA.fastq'
    barcode = 'AAAA' #in future extract this from filename
    tn_end = 'ACAGGTTG'
    trim_fastq(fastq_file, barcode, tn_end)

if __name__ == "__main__":
    main()
