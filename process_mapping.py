#!/usr/bin/env python
"""
Counts the bowtie hits at each position in each sample

Future - filter on the 16/17 bp positions before this step; maybe even before bowtie mapping
- change d to OrderedDict to keep contigs in order?

"""

import sys
import re
from collections import Counter

def multifasta2dict(fna):
    """
    Returns a dictionary of fasta headers and sequences

    Input file in multifasta format
    >header1
    agctag
    >header2
    gattta

    Function returns the fasta data as a dictionary
    {'header1': 'agctag', 'header2': 'gattta'}
    """
    d = {}
    with open(fna, 'r') as fi:
        header = ''
        first_line = False # next line is first line of sequence in contig
        for line in fi:
            line = line.rstrip()
            if line.startswith('>'): # header line; therefore new contig
                header = line[1:] # new header
                first_line = True
            else: # sequence line
                if first_line:
                    d[header] = line # create dictionary entry
                    first_line = False
                # would this work in one step?
                else:
                    appended_sequence = d[header] + line
                    d[header] = appended_sequence
        return d

def TA_sites(fna):
    """ Enumerate TA dinucleotides in a fasta nucleotide file

    Returns a dictionary of 5'-TA positions for each contig in a multifasta file.
    Assumes circular genome; i.e. will check if the last [-1] nucleotide is a T
    and the first [0] is an A

    Function calls multifasta2dict() and returns the nucleotide positions
    of TA sites in each contig as dictionary values:
    {'header1': [20, 22, 24, ... 1401104], 'header2': [5, 16, 24, ... 39910]}
    """
    di = multifasta2dict(fna)
    do = {}
    for header in di:
        sequence = di[header]
        do[header] = []
        i = 0   # nucleotide index
        # Checks the linear molecule (i.e. all except the last base)
        while i < len(sequence)-1:
            if sequence[i:i+2] == 'ta':
                do[header].append(i+1)
            i += 1
        # Checks last base circular around to first base
        if i == len(sequence)-1:
            if (sequence[-1] + sequence[0]) == 'ta':
                do[header].append(i+1)
    print(do)


def count_bowtie(bowtie_output):
    """ List counts of reads to the Left and Right of each TA site

    future:
    default show_all=False
    show_all=True will show all of the sites even those will count=0

    THIS MODULE STILL NEEDS TO BE WRITTEN!!!!!!!!

    """
    with open(bowtie_output, 'r') as fi:
        outfile = bowtie_output + '_counted'
        with open(outfile, 'w') as fo:
            count = {}
            for line in fi:
                sample_assignment = line.split('\t')[0].rfind('//')
                read_data_to_count = ('{0}\t{1}\t{2}\t{3}'.format(
                    line.split('\t')[0][sample_assignment+2:],
                    line.split('\t')[1],
                    line.split('\t')[2],
                    line.split('\t')[3]))
                count[read_data_to_count] = count.get(read_data_to_count, 0) + 1
            for c in count:
                print('{0}\t{1}'.format(c, count[c]))

def normalize_cpm(bowtie_output):
    """ insert """
    pass

def map_to_gene(normalized_output):
    """ insert """
    pass


# ===== Start here ===== #

def main():
    fna = sys.argv[1]
    TA_sites(fna)

if __name__ == "__main__":
    main()
