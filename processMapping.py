#!/usr/bin/env python
"""
Counts the bowtie hits at each position in each sample

Future - filter on the 16/17 bp positions before this step; maybe even before bowtie mapping
- change d to OrderedDict to keep contigs in order?

"""

import sys
import re
import csv
import collections
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

def insertionNucleotides(bowtieOutput, experiment=''):
    """ Define the TA dinucleotide of the insertion from the bowtie output

    Note bowtie output, where bnt = bowtie nucleotide
    + strand (sequence to L of Tn):
    <bnt><chromosomeSequence><TA> (TA is part of the sequence)
    - strand (sequence to R of Tn):
    <bnt><TA><chromosomeSequence> (TA is part of the sequence)

    Returns a list of tuples for the insertion nucleotides and orientations.
    Each insertion read is listed separately:
    experiment, barcode, contig, TAnucleotide, orientation
    [('Exp002', 'AAAA', 'contig1', 514141, 'R'), ('Exp003', 'GCTA', 'contig5', 8141, 'L')]

    Print output tab-delimited in a file called <experiment>_bowtieOutput_InsertionDetail.txt

    """
    with open(bowtieOutput, 'r') as fi:
        # Create new output document or overwrite previous document

        ## Pass the temp directory here in the future?
        ## Make this consistent in all of the scripts
        ## Also can only pass experiment name and recreate the rest.

        with open('{0}/temp/{0}_bowtieOutput_InsertionDetail.txt'.format(experiment), 'w') as fo:
            pass
        #insertions = []
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
                insertionNt = bnt + readLength - 1
                orient = 'L'
            else: # negative strand read
                insertionNt = bnt + 1
                orient = 'R'
            newInsertion = (experiment, sample, contig, insertionNt, orient)
            #insertions.append(newInsertion)
            with open('{0}/temp/{0}_bowtieOutput_InsertionDetail.txt'.format(experiment), 'a') as fo:
                fo.write('\t'.join(str(item) for item in newInsertion) + '\n')
    #return insertions

def insertionNucleotidesCount(experiment=''):
    """ List counts for each transposon insertion:

    NOT PRECISE ANYMORE -- UPDATE THIS DOCSTRING

    Returns a list of tuples in which the frequence of each orientation in the
    dataset is listed:
    experiment, barcode, contig, TAnucleotide, orientation, countL, countR, countTotal
    [('Exp002', 'AAAA', 'contig1', 514141, 20, 14, 34),
    ('Exp003', 'GCTA', 'contig5', 8141, 0, 1, 1)]

    Print output tab-delimited in a file called <experiment>_bowtieOutput_InsertionDetailCount.txt

    """

    #with open('{0}/temp/{0}_bowtieOutput_InsertionDetail.txt'.format(experiment)) as fi:
    #    fi_tsv = csv.reader(f, delimiter='\t')
    #    headers = next(fi_tsv)
    #    for row in fi_tsv:
    #        print(row)
    #
    #Currently giving IOError: [Errno 2] No such file or directory


    # Count the insertions (unique per experiment/sample/contig/nucleotide/orientation)

    """counter = collections.Counter()

    # read the raw insertion data.
    # Count identical experiment/sample/contig/nucleotide/orientation
    with open('{0}/temp/{0}_bowtieOutput_InsertionDetail.txt'.format(experiment), 'r') as fi:
        for line in fi:
            # create a tuple of the unique hit data
            tupi = tuple(line.rstrip().split('\t'))
            counter.update([tupi])

    # Write to tab-delimited file
    # experiment/sample/contig/nucleotide/orientation/count
    with open('{0}/temp/{0}_bowtieOutput_InsertionDetailCount__Temp.txt'.format(experiment), 'w') as fo:
        for tup in counter:
            for x in tup:
                fo.write('{0}\t'.format(x.strip()))
            fo.write('{0}\n'.format(counter[tup]))"""

    # read the raw count data: experiment/sample/contig/nucleotide/orientation/counts
    with open('{0}/temp/{0}_bowtieOutput_InsertionDetailCount__Temp.txt'.format(experiment), 'r') as fi2:
        # dictionary of counted tuples
        # {('Exp002', 'AAAA', 'contig1', 514141): (20, 14, 34)}
        #       Key = experiment/sample/contig/nucleotide
        #       Value = Lcounts/Rcounts/totalCounts
        do = {}
        for line in fi2:
            Lcounts = Rcounts = 0
            tupi = tuple(line.rstrip().split('\t'))
            tnOrientation = tupi[4]
            if tnOrientation == 'L':
                Lcounts = int(tupi[5])
            if tnOrientation == 'R':
                Rcounts = int(tupi[5])
            # e is a single entry of experiment/sample/contig/nucleotide
            e = tupi[0:4]
            if e in do:
                do[e] = (do[e][0]+Lcounts, do[e][1]+Rcounts, do[e][2]+Lcounts+Rcounts)
            else:
                do[e] = (Lcounts, Rcounts, Lcounts+Rcounts)

    # Write to tab-delimited file
    # experiment/sample/contig/nucleotide/Lcount/Rcount/totalCount
    with open('{0}/temp/{0}_bowtieOutput_InsertionDetailCount.txt'.format(experiment), 'w') as fo2:
        for e in do:
            for x in e:
                fo2.write('{0}\t'.format(x.strip()))
            for y in do[e]:
                fo2.write('{0}\t'.format(str(y)))
            fo2.write('\n')



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
