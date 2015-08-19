#!/usr/bin/env python
"""
Counts the bowtie hits at each position in each sample

Future - change d to OrderedDict to keep contigs in order?

"""

import sys
import re
import csv
import collections
from operator import itemgetter

# Not currently called in the pipeline.
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

# Not currently called in the pipeline.
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

    Prints a tab-delimited file {experiment}/temp/{experiment}_bowtieOutput_InsertionDetailCount.txt
    in which the count of each insertion, for each orientation (left/right) is listed
    as follows (unsorted):
    Exp002  AAAA  contig1  514141  20  14  34
    Exp003  GCTA  contig5  8141  0  1  1
    Columns represent:
    experiment | barcode | contig | TAnucleotide | countL | countR | countTotal

    ---

    As an intermediate step the data are first written to file
    {experiment}/temp/{experiment}_bowtieOutput_InsertionDetailCount__Temp.txt'
    in the following format (unsorted):
    Exp002  AAAA  contig1  514141  L  20
    Exp002  AAAA  contig1  514141  R  14

    Print output tab-delimited in a file called <experiment>_bowtieOutput_InsertionDetailCount.txt

    """

    # TODO: Use csv.reader for more robust tab-delimited handling
    #with open('{0}/temp/{0}_bowtieOutput_InsertionDetail.txt'.format(experiment)) as fi:
    #    fi_tsv = csv.reader(f, delimiter='\t')
    #    headers = next(fi_tsv)
    #    for row in fi_tsv:
    #        print(row)
    #
    #Currently giving IOError: [Errno 2] No such file or directory

    counter = collections.Counter()

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
            fo.write('{0}\n'.format(counter[tup]))

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

def filterSortCounts(experiment=''):
    """ Filter for min 3 reads total (1 each L/R) and maximum 10-fold L/R differential

    Also sort by the first four fields

    Data that do not pass this filter are treated as artefacts and not
    used in subsequent steps, such as data normalization.

    """
    with open('{0}/temp/{0}_bowtieOutput_InsertionDetailCount.txt'.format(experiment), 'r') as fi:
        li = []
        for line in fi:
            # experiment/sample/contig/nucleotide/Lcount/Rcount/totalCount
            tupi = tuple(line.rstrip().split('\t'))
            li.append(tupi)
        #print(li)

        # FILTER COUNT DATA
        lo = []
        for l in li:
            experiment = l[0]
            samples = l[1]
            contig = l[2]
            nucleotide = l[3]
            Lcounts = l[4]
            Rcounts = l[5]
            totalCounts = l[6]
            # minimum 3 total reads and minimum 1 read in each direction.
            if Lcounts >= 1 and Rcounts >= 1 and totalCounts >= 3:
                # maximum 10-fold L/R differential
                # L=1 R=10 ok, but not L=1 R=11
                if not (11 * min(Lcounts, Rcounts)) < (totalCounts):
                    lo.append(l)

        # SORT COUNT DATA
        loSorted = sorted(lo, key=itemgetter(0,1,2,3))

        # Write sorted/filtered data to tab-delimited file
        # experiment/sample/contig/nucleotide/Lcount/Rcount/totalCount
        with open('{0}/temp/{0}_bowtieOutput_InsertionDetailCountFilteredSorted.txt'.format(experiment), 'w') as fo:
            for e in loSorted:
                for x in e:
                    fo.write('{0}\t'.format(x.strip()))
                fo.write('\n')

def insertionCounts(experiment=''):
    """ Returns a list of tuples of insertion counts.

    [
    ('Exp002', 'AAAA', 'contig1', 514141, 20, 14, 34),
    ('Exp003', 'GCTA', 'contig2', 1441, 4, 5, 9)
    ]

    Elements represent:
    [0] experiment
    [1] barcode
    [2] contig
    [3] TAnucleotide
    [4] countL
    [5] countR
    [6] countTotal


    Currently this module reads in the tab-delimited CSV file generated by
    filterSortCounts() to serve as a basis for future analyses without
    the need to writing everything to a file (and then reading it back).

    A future, optimized version of the software might encasulate some of the
    previous manipulations directly in this module.

    """
    with open('{0}/temp/{0}_bowtieOutput_InsertionDetailCountFilteredSorted.txt'.format(experiment), 'r') as fi:
        counts = []
        for line in fi:
            tupi = tuple(line.rstrip().split('\t'))
            counts.append(tupi)
        return counts

def normalizeCpm(experiment):
    """ Normalize every sample to 1E6 CPM

    """
    # list of tuples of each mapped insertion, counted
    counts = insertionCounts(experiment)

    # Total count of filtered, mapped reads for each sample
    ddTotalCountBySample = collections.defaultdict(int)

    for tupi in counts:
        #Data for each insertion
        experiment = tupi[0]
        barcode = tupi[1]
        readCount = int(tupi[6])
        # Add to the total counts for that sample (barcode)
        ddTotalCountBySample[(experiment, barcode)] += readCount

    # Print counts for logging
    for entry in ddTotalCountBySample:
        experiment, barcode = entry
        print('{exp}\t{bc}\t{counts}'.
            format(exp=experiment, bc=barcode, counts=str(ddTotalCountBySample[entry])))

    # Transform the original counts data (total only) by dividing each count by
    # the total for the sample and multiplying by 1E6.
    # resulting dataset normCountsTotal. For each tuple:
    #[0] = experiment
    #[1] = barcode
    #[2] = contig
    #[3] = TAnucleotide
    #[4] = countTotal (from tupi[6] totalCounts)
    #  Note that Left and Right counts are not carried through.
    normCountsTotal = []
    for tupi in counts:
        rawCounts = int(tupi[6])
        sample = (tupi[0], tupi[1])
        totalSampleCounts = ddTotalCountBySample[sample]
        normCounts = float(1E6) * rawCounts / totalSampleCounts
        newTup = (tupi[0], tupi[1], tupi[2], tupi[3],
                    rawCounts, normCounts)
        print(newTup)

    return normCountsTotal



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
