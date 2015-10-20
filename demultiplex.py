#!/usr/bin/env python
"""
Demultiplexes a FASTQ file into multiple files by 5' end barcode

Output path includes the experiment and sample name:
(pysinseq)/samples/{experiment}/{sample}.fastq

"""

from utils import *
import re # I THINK I CAN DELETE THIS
import os # I THINK I CAN DELETE THIS
import gzip
import sys
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
    # For each line in the FASTQ file:
    #   Assign to barcode (fastq record into the dictionary)
    #   Then write to the appropriate output file
    with screed.open(fastq_file) as seqfile:
        for read in seqfile:
            #TODO: generalize for other length barcodes using
            # the length of the barcode and slice()
            try:
                barcode = read.sequence[0:4]
                demultiplex_dict[barcode].append(read)
            except:
                demultiplex_dict['_other'].append(read)
            # Every 10^6 sequences write and clear the dictionary
            nreads += 1
            if nreads % 1E6 == 0:
                writeReads(demultiplex_dict, barcodes_dict, experiment)
                # Clear the dictionary after writing to file
                for sampleName in demultiplex_dict:
                    demultiplex_dict[sampleName] = []
                sys.stdout.write('\r' + 'Records processed ... {:,}'.format(nreads))
    writeReads(demultiplex_dict, barcodes_dict, experiment)
    sys.stdout.write('\r' + 'Records processed ... {:,}'.format(nreads))

def writeReads(demultiplex_dict, barcodes_dict, experiment):
    """
    Write the fastq data to the correct (demultiplexed) file
    """
    for sampleName, barcode in barcodes_dict.items():
        with gzip.open('samples/{experiment}/{sampleName}.fastq.gz'.format( \
            experiment = experiment,
            sampleName = sampleName), 'a') as fo:
            for fastqRead in demultiplex_dict[barcode]:
                fo.write(bytes('@{n}\n{s}\n+\n{q}\n'.format( \
                    n = fastqRead.name,
                    s = fastqRead.sequence,
                    q = fastqRead.quality), 'UTF-8'))

def samplesToProcess(sample_file, experiment):
    """
    Returns a list of the sample paths to process in the current analysis.

    e.g.:
    ['samples/E001/E001_01.fastq.gz', 'samples/E001/E001_02.fastq.gz']

    """
    barcodes_dict = barcodes_prep(sample_file)
    sampleFile_list = []
    for sample in sorted(barcodes_dict):
        sampleFile_list.append('samples/{experiment}/{sample}.fastq.gz'.format(
            experiment=experiment,
            sample=sample
            ))
    return sampleFile_list

# ===== Start here ===== #

def main():
    fastq_file = sys.argv[1]
    sample_file = sys.argv[2]
    demultiplex_fastq(fastq_file, sample_file, 'E1001')

if __name__ == "__main__":
    main()
