#!/usr/bin/env python
"""
Demultiplexes a FASTQ file into multiple files by 5' end barcode.

Output path includes the experiment and sample name:
(pysinseq)/{experiment}/{sample}.fastq

"""

import collections
import csv
import gzip
import screed
import sys
from utils import convert_to_filename


def sample_prep(sample_file, barcode_qc):
    """
    Return ordered dictionary of sample name and barcode for each sample.

    samples = OrderedDict([('name1', {'name': 'name1', 'barcode': 'barcode1'}),
        ('name2', {'name': 'name2', 'barcode': 'barcode2'})])
    ignores comment lines in sample file that begin with #

    Exit if duplicate sample names or barcodes are identified.
    If barcode_qc=True, exit if duplicate barcodes are identified.
    """
    sampleDict = collections.OrderedDict()
    with open(sample_file, 'r', newline='') as csvfile:
        sampleReader = csv.reader(csvfile, delimiter='\t')
        for line in sampleReader:
            if not line[0].startswith('#'):
                # sample into a string that can be a filename; barcode to uppercase
                new_sample = convert_to_filename(line[0])
                if new_sample in sampleDict:
                    sys.stdout.write('Error: redundant sample identifier {}'.format(new_sample))
                    exit(1)
                try:
                    new_barcode = line[1].upper()
                except Exception:
                    # Of for downstream only if samples are already demultiplexed
                    new_barcode = ''
                else:
                    if barcode_qc:
                        if new_barcode == '':
                            sys.stdout.write('Missing barcode for sample {}'.format(new_sample))
                            exit(1)
                        if new_barcode in sampleDict.values():
                            sys.stdout.write('Error: redundant barcode {}'.format(new_barcode))
                            exit(1)
                sampleDict[new_sample] = {
                    'name': new_sample,
                    'barcode': new_barcode
                }
    return sampleDict


def demultiplex_fastq(fastq_file, samples_dict, experiment):
    """Demultiplex a fastq input file by 5' barcode into separate files."""
    barcodes_dict = {}
    for sample in samples_dict:
        barcodes_dict[sample] = samples_dict[sample]['barcode']
    # holder for unassigned barcodes
    barcodes_dict['_other'] = '_other'

    # Dictionary of lists to hold FASTQ reads until they are written to files
    demultiplex_dict = {}

    for sampleName, barcode in barcodes_dict.items():
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
            # TODO(generalize for other length barcodes using
            # the length of the barcode and slice()
            try:
                barcode = read.sequence[0:4]
                demultiplex_dict[barcode].append(read)
            except Exception:
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
    sys.stdout.write('\r' + 'Records processed ... {:,}'.format(nreads) + '\n')


def writeReads(demultiplex_dict, barcodes_dict, experiment):
    """Write the fastq data to the correct (demultiplexed) file."""
    for sampleName, barcode in barcodes_dict.items():
        with gzip.open('{experiment}/raw_data/{sampleName}.fastq.gz'.format(
            experiment=experiment,
            sampleName=sampleName), 'a'
        ) as fo:
            for fastqRead in demultiplex_dict[barcode]:
                fo.write(bytes('@{n}\n{s}\n+\n{q}\n'.format(
                    n=fastqRead.name,
                    s=fastqRead.sequence,
                    q=fastqRead.quality), 'UTF-8'))


def demultiplexedSamplesToProcess(sample_file, experiment):
    """
    Returns a list of tuples of the sample names and sample paths to process in the current analysis.

    e.g.:
    [('E001_01', 'experiment01/raw_data/E001_01.fastq.gz'),
        ('E001_02', 'experiment01/raw_data/E001_02.fastq.gz')]
    """
    barcodes_dict = sample_prep(sample_file, False)
    sample_list = []
    samplePath_list = []
    for sample in sorted(barcodes_dict):
        samplePath = '{experiment}/raw_data/{sample}.fastq.gz'.format(
            experiment=experiment,
            sample=sample
        )
        sample_list.append(sample)
        samplePath_list.append(samplePath)
    return sample_list, samplePath_list

"""
trim_fastq() and writeTrimmedReads() are derived from functions above.
They might be able to be consolidated into more general code. In utils.py?

"""


def trim_fastq(fastq_file, output_file, sampleName, bLen):
    """
    Write a fastq file with only the chromosome sequence.

    Writes a new fastq.gz file with:
    '_trimmed' appended to its name. 'sample01.fastq.gz' > 'sample01_trimmed.fastq.gz'
    (always done in gzipped format)
    Sequence/quality info for barcode and transposon data are not written
    into new file.
    bLen = barcode length at 5' end of read
    """
    # list to hold trimmed FASTQ reads until they are written to files
    trimmed_list = []
    # count of reads
    nreads = 0
    # For each line in the FASTQ file:
    #   Check for intact transposon.
    #   (Already know that barcode is ok if not in '_other' file)
    #   Then write to the appropriate output file
    with screed.open(fastq_file) as seqfile:
        for read in seqfile:
            try:
                # Identify length of chromosome sequence in read
                # (after barcode, before Tn)
                # Tn location as number of nuceotides after the barcode
                # Note that the 'TA' are in the chromosome too, so add 2
                chromosomeSeq = read.sequence.find('TAACAGGTTG') + 2 - bLen
                # slice of sequence and quality to extract
                seqSlice = slice(bLen, (chromosomeSeq + bLen))
                # Good read!
                if chromosomeSeq in range(16, 18):
                    trimmed_list.append(read[seqSlice])
            except Exception:
                pass
            # Every 10^6 sequences write and clear the dictionary
            nreads += 1
            if nreads % 1E6 == 0:
                writeTrimmedReads(trimmed_list, sampleName, output_file)
                # Clear the dictionary after writing to file
                trimmed_list = []
                sys.stdout.write('\r' + 'Records processed ... {:,}'.format(nreads))
    writeTrimmedReads(trimmed_list, sampleName, output_file)
    trimmed_list = []
    sys.stdout.write('\r' + 'Records processed ... {:,}'.format(nreads) + '\n')


def writeTrimmedReads(trimmed_fastq_list, sampleName, trimmed_fastq_filepath):
    """Write the trimmed fastq data into the experiment directory."""
    with open(trimmed_fastq_filepath, 'a') as fo:
        for fastqRead in trimmed_fastq_list:
            fo.write('@{n}\n{s}\n+\n{q}\n'.format(
                n=fastqRead.name,
                s=fastqRead.sequence,
                q=fastqRead.quality))


def main():
    """Start here."""
    # fastq_file = sys.argv[1]
    # sample_file = sys.argv[1]
    # demultiplex_fastq('_exampleData/example01.fastq',
    #    '_exampleData/example01.txt', 'example01')
    # trim_fastq('samples/example01/E001_01.fastq.gz')

if __name__ == "__main__":
    main()
