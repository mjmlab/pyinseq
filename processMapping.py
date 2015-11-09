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

def mapSites(bowtieOutput):
    pass

def trimBarcodeTransposon(fastq_path):
    """
    Trim barcode and transposon sequence, write to new file

    Assumes demultiplexed file. Barcode information is not retained.

    """

    # List to hold FASTQ reads until written to file
    trimmedReads_list = []
    # Clear the file to be written into
    clearTempReadFile()

    # Assume 4 bp 5' barcode for now
    bLen = 4

    # Length of chromosomal sequence
    # Should be 16-17 bp (filter on this)
    # Also used to match quality with corresponding sequence
    chromosomeSeq = 0

    # count of reads
    nreads = 0
    # For each line in the FASTQ file:
    #   Rewrite without barcode or transposon
    #   Then write to the appropriate output file
    with screed.open(fastq_path) as seqfile:
        for read in seqfile:
            # TODO: Improve speed by writing own parser that slices out
            #    the barcode
            #    https://screed.readthedocs.org/en/v0.9/screed.html#writing-custom-sequence-parsers
            # intact transposon junction with 16-17 bp chromosomal DNA between
            # barcode and transposon

            # Tn location as number of nuceotides after the barcode
            # Note that the 'TA' are in the chromosome too, so add 2
            chromosomeSeq = read.sequence.find('TAACAGGTTG') + 2 - bLen

            #slice of sequence and quality to extract
            seqSlice = slice(bLen:(chromosomeSeq+bLen))

            # Good read!
            if chromosomeSeq in range(16, 18):
                # screed will cleverly use the sliced part of the sequence and
                # quality but the entire name/identifier
                trimmedReads_list.append(read[seqSlice])
            # ignore if does not pass filters above
            else:
                pass

            # Every 10^6 sequences write and clear the list
            nreads += 1
            if nreads % 1E6 == 0:
                writeTempReads(trimmedReads_list)
                # Clear the list after writing to file
                trimmedReads_list = []
                sys.stdout.write('\r' + 'Records processed ... {:,}'.format(nreads))
    writeReads(demultiplex_dict, barcodes_dict, experiment)
    sys.stdout.write('\r' + 'Records processed ... {:,}'.format(nreads))

def clearTempReadFile():

    # TODO: experiment should be a global variable. How to do that?

    """
    Clear the temporary fastq file that is used for bowtie input.
    """
    try:
        with open('{experiment}/temp/temp_bowtie_file.fastq'.format(
            experiment='E001'
            ), 'w') as fo:
            pass
    except FileNotFoundError:
        'No temp directory in which to create the bowtie input file.'

def writeTempReads(fastq_list):
    """
    Write the fastq data to the correct (demultiplexed) file
    """
    for sampleName, barcode in barcodes_dict.items():
        with open('temp/temp_bowtie_file.fastq', 'a') as fo:
            for fastqRead in demultiplex_dict[barcode]:
                fo.write(bytes('@{n}\n{s}\n+\n{q}\n'.format( \
                    n = fastqRead.name,
                    s = fastqRead.sequence,
                    q = fastqRead.quality), 'UTF-8'))


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
    """ Filter for min 1 L read, 1 R read and maximum 10-fold L/R differential

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
            experiment, samples, contig, nucleotide, Lcounts, Rcounts, \
                totalCounts = l[0:7]
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

    Returns a list of tuples:

    [
    ('Exp001', 'TTTT', 'contig1', '999401', 80.00268809031984)
    ]

    experiment, barcode, contig, TAnucleotide, normalized_counts
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
    sample_file = sys.argv[2]
    experiment = sys.argv[3]
    processList = samplesToProcess(sample_file, experiment)
    for s in processList:
        print(s)

if __name__ == '__main__':
    main()
