#!/usr/bin/env python
"""
Counts the bowtie hits at each position in each sample

Future - filter on the 16/17 bp positions before this step; maybe even before bowtie mapping
- change d to OrderedDict to keep contigs in order?

"""

import sys
import re
from collections import Counter
from operator import itemgetter

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
        firstLine = False # next line is first line of sequence in contig
        for line in fi:
            line = line.rstrip()
            if line.startswith('>'): # header line; therefore new contig
                header = line[1:] # new header
                firstLine = True
            else: # sequence line
                if firstLine:
                    d[header] = line # create dictionary entry
                    firstLine = False
                # would this work in one step?
                else:
                    appendedSequence = d[header] + line
                    d[header] = appendedSequence
        return d

def taSites(fna):
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

def insertionNucleotides(bowtieOutput):
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
    with open(bowtieOutput, 'r') as fi:
        insertions = []
        for line in fi:
            # Redo this with regex
            idLine = line.rstrip().split('\t')[0]
            loc1 = idLine.rfind('//')
            loc2 = idLine[:idLine.rfind(':')].rfind(':')
            loc3 = idLine.rfind(':')
            experiment = idLine[loc1+2:loc2]
            sample = idLine[loc2+1:loc3]
            contig = line.rstrip().split('\t')[2]
            # Read direction
            if line.rstrip().split('\t')[1] == '+':
                direction = True # positive strand read
            else:
                direction = False # minus strand read
            # Calculate transposon insertion point
            readLength = len(line.rstrip().split('\t')[4])
            bnt = int(line.rstrip().split('\t')[3]) # bowtie nucleotide
            if direction: # positive strand read
                insertionNt = bnt + readLength - 1 # -2 if I made a +1 error
                orient = 'L'
            else: # negative strand read
                insertionNt = bnt + 1 # delete +1 if I made a +1 error
                orient = 'R'
            newInsertion = (experiment, sample, contig, insertionNt, orient)
            insertions.append(newInsertion)
    return insertions

def insertionCounts(bowtieOutput, experiment=''):
    """ List counts for each transposon insertion:

    Returns a list of tuples in which the frequence of each orientation in the
    dataset is listed:
    experiment, barcode, contig, TAnucleotide, orientation, countL, countR, countTotal
    [('Exp002', 'AAAA', 'contig1', 514141, 20, 14, 34),
    ('Exp003', 'GCTA', 'contig5', 8141, 0, 1, 1)]

    Print the output to a file called <experiment>_bowtieOutputSummarized.txt

    """
    # list - each insertion read listed individually
    li = insertionNucleotides(bowtieOutput)

    # list of hits in bowtie results independent of orientation
    # add placeholder integers for L, R, and Total counts
    uniqueHits = list(set([tuple(x for x in y[0:4]) for y in li]))

    # dictionary of counts of L and R counts
    # key is the index in the uniqueHits list of tuples
    counts = {}
    for i in range(0,len(uniqueHits)):
        counts[i] = [0, 0, 0]

    # loop through the bowtie output of the insertions
    # for each insertion...
    # check the index of the read in the uniqueHits list of tuples
    # and incremement the L (or R) and Total counts
    for i in li:
        hitsIndex = uniqueHits.index(i[0:4])
        if i[4] == 'L':
            counts[hitsIndex][0] += 1  # L
            counts[hitsIndex][2] += 1  # Total
        if i[4] == 'R':
            counts[hitsIndex][1] += 1  # R
            counts[hitsIndex][2] += 1  # Total

    lo = []
    for i in range(0,len(uniqueHits)):
        lo.append(uniqueHits[i] + tuple(counts[i]))

    with open('{0}_bowtieOutputSummarized.fastq'.format(experiment), 'w') as fo:
        # Sort list by first 4 columns
        for i in sorted(lo, key=itemgetter(0,1,2,3)):
            fo.write('\t'.join(str(item) for item in i) + '\n')
    #return(lo)

def filterCounts(bowtieOutput):
    """ Filter for min 5 reads total (1 each L/R) and maximum 10-fold L/R differential

    """
    li = insertionCounts(bowtieOutput)
    lo = []
    for i in li:
        # minimum 5 total reads and minimum 1 read in each direction.
        if i[4] >= 1 and i[5] >= 1 and i[6] >= 5:
            # maximum 10-fold L/R differential
            # L=1 R=10 ok, but not L=1 R=11
            if not (11 * min(i[4],i[5])) < (i[6]):
                print(i)
                lo.append(i)
    print(len(lo))

def normalizeCpm(bowtieOutput):
    """ Normalize every sample to 1E6 CPM

    START HERE!!

    """
    # list - each insertion read listed individually
    li = insertionCounts(bowtieOutput)
    totalCounts = [sum(i[6]) for i in zip(*li)]

def mapToGene(normalizedOutput):
    """ insert """
    pass


# ===== Start here ===== #

def main():
    bowtieOutput = sys.argv[1]
    experiment = sys.argv[2]
    insertionCounts(bowtieOutput, experiment)

if __name__ == '__main__':
    main()
