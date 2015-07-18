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

    Returns a list of tuples for the insertion nucleotides and orientations.
    Each insertion read is listed separately:
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
    return insertions

def insertion_counts(bowtie_output):
    """ List counts for each transposon insertion:

    Returns a list of tuples in which the frequence of each orientation in the
    dataset is listed:
    experiment, barcode, contig, TAnucleotide, orientation, countL, countR, countTotal
    [('Exp002', 'AAAA', 'contig1', 514141, 20, 14, 34),
    ('Exp003', 'GCTA', 'contig5', 8141, 0, 1, 1)]

    """
    # list - each insertion read listed individually
    li = insertion_nts(bowtie_output)
    # counted list - each insertion once
    lo = []

    #
    # Might be better to loop from nt 1 through << length of the genome >>
    # look for matches for left or right and count them up.
    # Maybe when making fna files store the lengths of the contigs then?
    #

    # Count the occurrance of each L or R read per contig/nt for that experiment/sample
    counts = {}
    for tup in set(li):
        counts[tup] = li.count(tup)

    for key in counts:
        print(key)
        print(counts[key])


    #s_li = sorted(counts)
    #for i in counts:
            #print(i)



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
    insertion_counts(fna)

if __name__ == '__main__':
    main()
