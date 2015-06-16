#!/usr/bin/python
# Process FASTQ Files for INSEQ
# Author: Mark J. Mandel
# Date: 2015-05-24


#class fastq()
#    def __init_(self, file):

import timing   # how long the program takes to run
import gzip
from itertools import islice
import os
import argparse


#input_file = ''
#barcode_file = ''

parser = argparse.ArgumentParser()
parser.add_argument('-i', help="input file")
args = parser.parse_args()
print args.echo

#    description='Processes INSeq data\n\n \
#    Currently the software uses barcodes in the barcodes dictionary \
#    and demultiplexes the input fastq file into separate files for each barcode.')
#parser.add_argument('-i', '--input', action='store', metavar=input_file,
#    required=True,
#    help='Input file in fastq or fastq.gz format.')
#parser.add_argument('-b', '--barcodes', action='store', destination=barcode_file,
#    default='barcodes.txt', required=True,
#    help='Tab-delimited list of barcodes and samples. See format in barcode_example.txt')


#   help
#   -i input
#   -s sample list
#   -g genome (fasta? ptt etc.)
#   -t transposon sequence OR...
#   -l transposon sequence left
#   -r transposon sequence right


class RawData(object):
    """RawData = the unprocessed FASTQ file and experimental info"""
    def __init__(self, fastq_file, sample_list):
        self.filepath = filepath
        self.sample_list = sample_list

class Dataset(object):
    """Dataset = each Barcoded Sample"""
    def __init__(self, filepath, barcode):
        self.barcode = barcode
        self.filepath = filepath

## MOVE THIS MODULE TO WITHIN RAWDATA ##




def bowtie_setup(bowtie_path, fna, ptt):
    """ Initializes the bowtie indexes using the provided .fna and .ptt files

    Concatenates the FASTA nucleotide (.fna) files into a single file
    in the case of multiple replicons in the genome, names the FASTA headers
    to correspond to the filenames in the protein features (.ptt) files, and
    then indexes the genome using bowtie-build"""
    pass


def demultiplex_fastq(fastq_file):
    """ Demultiplexes a FASTQ file into multiple files by 5' end barcode

    XXXX Write here about how the files get named and written.
    Barcodes remain in file at this point so that these files can be
    deposited in data repositories (e.g., GenBank).
    Discards sequences that not have a perfect match to the tranpsoson.
    XXXX
    """
    print('Parsing FASTQ input file')
    if fastq_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    barcodes = { \
    'GTAC':'Input1','AAAA':'Input2', \
    'TATA':'Input3','GTCA':'Input4', \
    'TTAA':'Output5','AACC':'Output6', \
    'GCTA':'Output7','AGTC':'Output8'}

    with opener(fastq_file, 'r') as f:

        # load barcode file


        "barcode length"
        b_len = 4

        """fastq record
        record[0] #identifier
        record[1] #sequence
        record[2] #alt_identifier
        record[3] #quality
        """
        record = []

        for i,line in enumerate(f):
            record.append(line.rstrip('\n'))
            if i % 4 == 3:
                if record[1][0:b_len] in barcodes:
                    with open('__{0}.fastq'.format(record[1][0:b_len]), 'a') as fo:
                        fo.write('\n'.join(record) + '\n')
                record = []

                """ Add a counter to show how many records have been processed """


    """Call a separate module that creates the barcode files (blank)
    Then append to those files here. Read barcodes in as a dictionary and
    create files into which fastq data for that barcode will be written"""
    ### NEED TO REALLY READ IN BARCODES HERE ###


"""DO THIS IN OUTPUT SUBDIRECTORY?"""
"""NEXT STEPS == CREATE FILE ON FIRST INSTANCE, THEN APPEND ON FUTURE INSTANCES?"""
"""MOVE MUTANT BARCODE FILES TO SEPARATE DIRECTORY"""
"""SHOW NUMBER OF FILES PROCESSED IN COMMAND LINE"""



def filter_sequences():
    """ Removes barcode and transposon data from sequence file.

    Discards files without perfect matches to the transposon.
    (Sequences not matching a barcode are removed too, but in reality
    these are likley removed at demultiplexing.)"""
    pass

def map_sequences():
    """ Maps sequences with bowtie

    """
    pass

def map_to_genes():
    """ Maps hits to genes

    """
    pass

def compute_statistics():
    """ Computes log ratios and significance for change from input to output

    """
    pass

def list_results():
    """ Lists gene-level results in a tab-delimited spreadsheet

    Prints the log file to a text file called INSEQ_Experiment_Log.txt
    """

def plot_results():
    """ Plots key features of the dataset.

    """
    pass


demultiplex_fastq('testfile.fastq.gz')

""" Write statistics into a log file"""




#Experiment_E0001 = RawData('sample1.fastq.gz','E0001_sample_list.txt')

# in the __init__.py define the transposon sequence.

# test if lines 1, 5, 9 start with @ and if lines 2, 6, 10 have DNA;
#       otherwise exit and describe likely not the right format.
# for i in [0,4,8]
# Need to do the processing of the barcode and transposon
# Print some summary data - number of reads, number of reads per barcode
# Put the summary data in a separate output file
