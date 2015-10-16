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

from utils import *
import re
import gzip
import sys  # temporary - for command line arg
import os
import argparse
import screed

def barcodes_prep(sample_file):
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
    barcode_dict = {}
    with open(sample_file, 'r') as fi:
        for line in fi:
            if not line.startswith('#'):
                new_sample = line.rstrip().split('\t')[0]
                new_sample = convert_to_filename(new_sample)
                if new_sample in barcode_dict:
                    print('Error: redundant sample identifier {}'.format(new_sample))
                    exit(1)
                new_barcode = line.rstrip().split('\t')[1].upper()
                if new_barcode in barcode_dict.values():
                    print('Error: redundant barcode {}'.format(new_barcode))
                    exit(1)
                barcode_dict[new_sample] = new_barcode
                barcode_length = len(new_barcode)

    # Print barcodes and check for same length
    for b in sorted(barcode_dict):
        print('{0}\t{1}'.format(b, barcode_dict[b]))
        if len(barcode_dict[b]) != barcode_length:
            print('Error: barcodes are not the same length')
            exit(1)

    print('n={0} unique barcodes of length {1} nt'.format(len(barcode_dict), barcode_length))

    return barcode_dict

def demultiplex_fastq(fastq_file, sample_file, experiment):
    """
    Demultiplex a fastq input file by 5' barcode into separate files

    """

    #if fastq_file.endswith('.gz'):
    #    opener = gzip.open
    #else:
    #    opener = open

    barcodes_dict = barcodes_prep(sample_file)
    # holder for unassigned barcodes
    barcodes_dict['_other'] = '_other'

    # Dictionary of lists to hold FASTQ reads until they are written to files
    demultiplex_dict = {}

    for sampleName, barcode in barcodes_dict.items():
        # create the /pyinseq/samples/{experiment}/{sample} files
        # exit if any of the files already exist
        createDemultiplexFiles(experiment, sampleName)
        # create a dictionary to hold fastq data before it is written
        # key = barcode (NOT THE SAMPLE)
        # values = list of fastq records
        demultiplex_dict[barcode] = []

    # count of reads
    nreads = 0
    # Read in the FASTQ file
    # Assign to barcode
    # Write to the appropriate output file
    with screed.open(fastq_file) as seqfile:
        for read in seqfile:

            #TODO: generalize for other length barcodes
            try:
                barcode = read.sequence[0:4]
                demultiplex_dict[barcode].append(read)
            except:
                demultiplex_dict['_other'].append(read)
            # Every 10^5 sequences write and clear the dictionary
            nreads += 1
            if nreads % 10000 == 0:
                writeReads(demultiplex_dict, barcodes_dict, experiment)
                print(nreads, 'records processed')
    writeReads(demultiplex_dict, barcodes_dict, experiment)
    print(nreads, 'records processed')

def writeReads(demultiplex_dict, barcodes_dict, experiment):
    """
    Write the fastq data to the correct (demultiplexed) file
    """
    for sampleName, barcode in barcodes_dict.items():
        with open('samples/{experiment}/{sampleName}.fastq'.format( \
            experiment = experiment,
            sampleName = sampleName), 'a') as fo:
            for fastqRead in demultiplex_dict[barcode]:
                fo.write('@{n}\n{s}\n+\n{q}\n'.format( \
                    n = fastqRead.name,
                    s = fastqRead.sequence,
                    q = fastqRead.quality))

#fo.write('@{n}\n{s}\n+{a}\n{q}\n'.format \

# ===== Start here ===== #

def main():
    fastq_file = sys.argv[1]
    sample_file = sys.argv[2]
    demultiplex_fastq(fastq_file, sample_file, 'E001')

if __name__ == "__main__":
    main()
