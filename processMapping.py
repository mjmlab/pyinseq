#!/usr/bin/env python
"""
Counts the bowtie hits at each position in each sample

"""

from demultiplex import barcodes_prep
import os
import sys
import re
import csv
import collections
import screed
from operator import itemgetter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def mapSites(bowtieOutput):
    # Placeholder for dictionary of mapped reads in format:
    # {(contig, position) : [Lcount, Rcount]}
    mapDict = {}
    with open(bowtieOutput, 'r') as fi:
        for line in fi:
            bowtiedata = line.rstrip().split('\t')
            # Calculate transposon insertion point = transposonNT
            readLength = len(bowtiedata[4])
            contig = str(bowtiedata[2])
            insertionNT = int(bowtiedata[3])
            if bowtiedata[1] == '+': # positive strand read
                insertionNt = insertionNT + readLength - 1
                mapDict.setdefault((contig,insertionNt),[0,0])[0] += 1   # Lcount
            else: # negative strand read
                insertionNt = insertionNT + 1
                mapDict.setdefault((contig,insertionNt),[0,0])[1] += 1   # Rcount
    # write tab-delimited of contig/nucleotide/Lcount/Rcount/TotalCount
    root, ext = os.path.splitext(bowtieOutput)
    with open('{0}_mapped{1}'.format(root, ext), 'a') as fo:
        writer = csv.writer(fo, delimiter='\t', dialect='excel')
        header_entry = ('contig', 'nucleotide', 'left_counts', 'right_counts', 'total_counts')
        writer.writerow(header_entry)
        for insertion in sorted(mapDict):
            row_entry = (insertion[0], insertion[1], mapDict[insertion][0], mapDict[insertion][1], mapDict[insertion][0] + mapDict[insertion][1])
            writer.writerow(row_entry)
    return(mapDict)

def mapGenes():
    pass


##  Not filtering now, showing all results
##  Will re-implement filtering once re-write these methods (e.g., in pandas)
def filterSortCounts(experiment, sample):
    """ Filter for min 1 L read, 1 R read and maximum 10-fold L/R differential

    Also sort by the first four fields

    Data that do not pass this filter are treated as artefacts and not
    used in subsequent steps, such as data normalization.

    """
    with open('{0}/{1}_output_bowtie_mapped.txt'.format(experiment, sample), 'r') as fi:
        li = []
        for line in fi:
            # contig/nucleotide/Lcount/Rcount/TotalCount
            tupi = tuple(line.rstrip().split('\t'))
            li.append(tupi)
        #print(li)

        # FILTER COUNT DATA
        lo = []
        for l in li:
            contig, nucleotide, Lcounts, Rcounts, totalCounts = l[0:5]
            lo.append(l)
            # TODO: RETURN TO FILTERING BELOW
            """# minimum 1 read in each direction.
            if int(Lcounts) >= 1 and int(Rcounts) >= 1:
                # maximum 10-fold L/R differential
                # L=1 R=10 ok, but not L=1 R=11
                if not (11 * min(Lcounts, Rcounts)) < (totalCounts):
                    lo.append(l)"""

        # SORT COUNT DATA
        loSorted = sorted(lo, key=itemgetter(0,1,2,3))

        # Write sorted/filtered data to tab-delimited file
        # experiment/sample/contig/nucleotide/Lcount/Rcount/totalCount
        with open('{0}/{1}_output_bowtie_mapped_filtered.txt'.format(experiment, sample), 'w') as fo:
            for e in loSorted:
                for x in e:
                    fo.write('{0}\t'.format(x.strip()))
                fo.write('\n')

def normalizeCpm(experiment):
    """ Normalize every sample to 1E6 CPM

    Returns a list of tuples:

    [
    ('contig1', '999401', 80.00268809031984)
    ]

    contig, nucleotide, normalized_counts
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

    # TODO: Write counts for logging
    """
    for entry in ddTotalCountBySample:
        experiment, barcode = entry
        print('{exp}\t{bc}\t{counts}'.
            format(exp=experiment, bc=barcode, counts=str(ddTotalCountBySample[entry])))
    """

    # Transform the original counts data (total only) by dividing each count by
    # the total for the sample and multiplying by 1E6.
    # resulting dataset normCountsAll. For each tuple:
    #[0] = experiment
    #[1] = barcode
    #[2] = contig
    #[3] = TAnucleotide
    #[4] = countTotal (from tupi[6] totalCounts)
    #  Note that Left and Right counts are not carried through.
    normCountsAll = []
    for tupi in counts:
        # note: can add back rawCounts if it would be valuable
        # to have them in the output
        rawCounts = int(tupi[6])
        sample = (tupi[0], tupi[1])
        totalSampleCounts = ddTotalCountBySample[sample]
        normCounts = float(1E6) * rawCounts / totalSampleCounts
        newTup = (tupi[0], tupi[1], tupi[2], tupi[3], normCounts)
        normCountsAll.append(newTup)
    return normCountsAll

def fttLookupTable(organism, experiment=''):
    """
    Import the ftt file and process as a lookup table

    No headers
    Includes the contig in every row
    Separates out the start..end to separate start, end

    """
    fttLookup = []
    with open('{0}/temp/{1}.ftt'.format(experiment, organism), 'r') as ftt:
        for line in ftt:
            # Capture contig information
            if line.startswith('LOCUS'):
                contig = line.rstrip().split('\t')[1]
            # Only print if first digit is numeric.
            # Note that on ptt files sometimes the first digit will be a
            # less-than sign (<) but currently we exclude those characters
            # in generating the .ftt file. Return to this criterion if necessary.
            # Alternate way would be first character is not 'L' (from Locus or Location)
            if unicode(line[0], 'utf-8').isnumeric():
                fttImport = line.rstrip().split('\t')
                # Split the location into separate elements:
                # 1014..4617 > 1014, 4617
                start, end = fttImport[0].split('..')[0:2]
                fttTemp = [contig] + [start] + [end] + fttImport[1:]
                fttLookup.append(fttTemp)
    return fttLookup

def mapToGene(organism, experiment=''):
    """
    Given a set of insertions with nucleotide data, provide the corresponding
    detail from the .ftt file:
    contig
    locus tag
    gene
    product
    beginning nucleotide
    ending nucleotide
    insertion location in gene (5' end = 0.0, 3' end = 1.0)

    Note that insertions that hit multiple (overlapping) genes will have
    multiple output lines printed.

    """

    fttLookup = fttLookupTable(organism, experiment)

    # Import the insertion data
    normCountsAll = normalizeCpm(experiment)
    mappedHitList = []
    for sample in normCountsAll:
        insertionContig = sample[2]
        insertionNucleotide = int(sample[3])
        # Used to save previous feature for intergenic calling
        prevFeature = ''
        for i, feature in enumerate(fttLookup):
            genomeContig = ''
            featureStart = 0
            featureEnd = 0
            genomeContig = feature[0]
            # number of genes/features hit.
            # 0 = intergenic. 1 = a gene. 2 = overlapping genes.
            # if 2, will print 2 entries in the output.
            genesHit = 0
            if insertionContig == genomeContig:
                featureStart, featureEnd = int(feature[1]), int(feature[2])
                orientation = feature[3]
                TAnucleotide = float(sample[3])
                if insertionNucleotide >= featureStart:
                    if insertionNucleotide <= featureEnd:
                        # 0.0 = 5'end ; 1.0 = 3'end

                        # TODO: Should featureEnd had +1 added?
                        if orientation == '+':
                            threePrimeness = (TAnucleotide - featureStart) / \
                                (featureEnd - featureStart)
                        if orientation == '-':
                            threePrimeness = (featureEnd - TAnucleotide) / \
                                (featureEnd - featureStart)
                        sampleDetails = [sample[0], sample[1], sample[2], \
                            sample[3], sample[4], threePrimeness]
                        featureDetails = [feature[6], feature[7], \
                            feature[10], feature[1], feature[2], \
                            feature[3], feature[4], feature[5]]

                        # experiment, barcode, contig, TAnucleotide, CPM,
                        # threePrimeness, gene, locusTag, description,
                        # featureStart, featureEnd, strand, length, PID
                        mappedHit = sampleDetails + featureDetails
                        mappedHitList.append(mappedHit)
                        # counts for overlap. Not being used currently
                        genesHit += 1

                    # TODO: Intergenic regions
                    # The logic here is only a start.
                    # Also note that intergenic hits at the end of a contig
                    # should be listed as between the last gene of the contig
                    # and the first gene of the *same* contig.
                    else:
                        pass
                        #if genesHit == 0
                        #    try:
                        #        print(fttLookup[i+1])  # the feature after (feature)
                        #    except:
                        #        print(fttLookup[0])
            prevFeature = feature
        genesHit = 0
    with open('{0}/temp/{0}_bowtieOutput_MappedToGeneDetail.txt'.format(experiment), 'w') as fo:
        for hit in mappedHitList:
            for entry in hit:
                fo.write('{0}\t'.format(entry))
            fo.write('\n')
    return mappedHitList

def mapToGeneSummary(cutoff, organism, experiment=''):
    """
    For each entry in a feature table (.ftt) list the summary of hits
    for each sample in the experiment

    Count only those hits that fall at or below the cutoff specified.
    e.g., cutoff=0.9 means that genes with threePrimeness <= 0.9 will be counted
    and hits in the 3-prime most ~10% of the gene will not be included.
    """

    featureTable = fttLookupTable(organism, experiment)

    # Header will be extended in the future to
    # list the experiment and barcode of each sample of interest
    baseHeader = ['Contig', 'Start', 'End', 'Strand', 'Length', 'PID',
        'Gene', 'Synonym', 'Code', 'COG', 'Product']

    # Add header row to ftt file
    featureTable.insert(0, baseHeader)

    # Get individual mapped hits
    mappedHitList = mapToGene(organism, experiment)

    ## IF barcode is in header[:-1] then check gene and threePrimeness.
    ## IF barcode is not in header, add it + [barcode]. Then check gene....
    ## For each gene do the math in a variable. Then add the result of that variable
    ## before moving on to the next feature. Every feature needs a 0 or a higher value.

    # - featureTable[0][-1] barcode.         For feature in table....:
    #Set up table with index, start, end.
    # Match for all in a barcode then add to entries. Repeat for next barcode.

    # current column, sample being matched
    # sample format is experiment:barcode
    currentColumn = len(featureTable[0])-1
    currentSample = ''

    for hit in mappedHitList:
        currentSample = featureTable[0][-1]
        sample = ('{exp}:{sample}'.format(exp=str(hit[0]), sample=str(hit[1])))
        threePrimeness = int(hit[5])
        if currentSample != sample:
            # Add the new sample (experiment:barcode) as a new column the table
            # and fill in entire column with 0's
            currentColumn += 1
            featureTable[0].append(sample)
            # Fill the rest of the new column with 0 as the count for each gene
            for feature in featureTable[1:]:
                feature.append(0)
        # Go through each hit in the mappedHitList for that barcode and
        # add the normalized counts to the respective gene if the hit falls
        # within the 5' cutoff
        if threePrimeness <= cutoff:
            for i, f in enumerate(featureTable):
                hitLocusTag = hit[7]
                fttLocusTag = f[7]
                # matches based on locusTag.
                # In future should I instead create an index field in the .ftt?
                if hit[7] == f[7]:
                    featureTable[i][currentColumn] += hit[4]
    with open('{0}/temp/{0}_bowtieOutput_MappedToGeneSummary.txt'.format(experiment), 'w') as fo:
        for line in featureTable:
            for entry in line:
                fo.write('{0}\t'.format(entry))
            fo.write('\n')


# ===== Start here ===== #

def main():
    bowtieOutput = sys.argv[1]
    mapSites(bowtieOutput)
    #sample_file = sys.argv[2]
    #experiment = sys.argv[3]
    #processList = samplesToProcess(sample_file, experiment)
    #for s in processList:
    #    print(s)

if __name__ == '__main__':
    main()
