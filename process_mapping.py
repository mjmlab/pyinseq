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
    return do


def insertion_nts(bowtie_output):
    """ Define the TA dinucleotide of the insertion from the bowtie output

    Note bowtie output, where bnt = bowtie nucleotide
    + strand (sequence to L of Tn):
    bnt_sequence_TA
    - strand (sequence to R of Tn):
    bnt_TA_sequence

    Returns a list of tuples for the insertion nucleotides and orientations:
    experiment, barcode, contig, TAnucleotide, orientation
    [('Exp002', 'AAAA', 'contig1', 514141, 'R'), ('Exp003', 'GCTA', 'contig5', 8141, 'L')]

    CHECK -- DID MY TEST FNA FILE SHIFT BY 1 NT?

    """
    with open(bowtie_output, 'r') as fi:
        insertions = []
        for line in fi:
            # Redo this with regex
            idline = line.rstrip().split('\t')[0]
            loc1 = idline.rfind('//')
            loc2 = idline[:idline.rfind(':')].rfind(':')
            loc3 = idline.rfind(':')
            experiment = idline[loc1+2:loc2]
            sample = idline[loc2+1:loc3]
            contig = line.rstrip().split('\t')[2]
            # Read direction
            if line.rstrip().split('\t')[1] == '+':
                direction = True # positive strand read
            else:
                direction = False # minus strand read
            # Calculate transposon insertion point
            read_length = len(line.rstrip().split('\t')[4])
            bnt = int(line.rstrip().split('\t')[3]) # bowtie nucleotide
            if direction: # positive strand read
                insertion_nt = bnt + read_length - 1 # -2 if I made a +1 error
                orient = 'L'
            else: # negative strand read
                insertion_nt = bnt + 1 # delete +1 if I made a +1 error
                orient = 'R'
            new_insertion = (experiment, sample, contig, insertion_nt, orient)
            insertions.append(new_insertion)
    print(insertions)
    return insertions

def count_bowtie(bowtie_output, fna):
    """ List counts for each TA site:

    {'header': [[site1, left_count, right_count, total_count], [site2, left_count, right_count, total_count]]}

    future:
    default show_all=False
    show_all=True will show all of the sites even those will count=0

    """
    ta_di = TA_sites(fna)   # dictionary of ta sites in each contig
    do = {}

    # For each TA site in ta_di {'header1': [site1, site2]} ...
    # Set up nested list in do {'header1': [[site1, 0, 0, 0], [site2, 0, 0, 0]]}



    """for header in ta_di:
        for ta

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
                print('{0}\t{1}'.format(c, count[c]))"""

def normalize_cpm(bowtie_output):
    """ insert """
    pass

def map_to_gene(normalized_output):
    """ insert """
    pass


# ===== Start here ===== #

def main():
    fna = sys.argv[1]
    bowtie_insertions(fna)

if __name__ == "__main__":
    main()
