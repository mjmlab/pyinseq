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
    # overallTotal = denominator for cpm calculation
    overallTotal = 0
    cpm = 0
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
            overallTotal += 1
    # write tab-delimited of contig/nucleotide/Lcount/Rcount/TotalCount/cpm
    # use the index totalSampleCounts as the denominator for cpm calculation
    root, ext = os.path.splitext(bowtieOutput)
    with open('{0}_mapped{1}'.format(root, ext), 'a') as fo:
        writer = csv.writer(fo, delimiter='\t', dialect='excel')
        header_entry = ('contig', 'nucleotide', 'left_counts', 'right_counts', 'total_counts', 'cpm')
        writer.writerow(header_entry)
        for insertion in sorted(mapDict):
            Lcounts = mapDict[insertion][0]
            Rcounts = mapDict[insertion][1]
            totalCounts = mapDict[insertion][0] + mapDict[insertion][1]
            cpm = float(1E6) * totalCounts / overallTotal
            row_entry = (insertion[0], insertion[1], Lcounts, Rcounts, totalCounts, cpm)
            writer.writerow(row_entry)
    return(mapDict)

def mapGenes(organism, sample, experiment=''):
    """
    Maps insertions to genes

    1. Writes a csv file listing the gene for each insertion with the tabs:
    contig, nucleotide, Lcounts, Rcounts, totalCounts, cpm, threePrimeness, locus_tag
    2. Returns a dictionary of aggregate counts per gene (for any gene with at
    least one hit)

    Insertions that map to multiple genes result in multiple lines in the CSV file
    and are counted for both genes in the returned dictionary.

    ThreePrimeness = insertion location in gene (5' end = 0.0, 3' end = 1.0)
    At present this is not used. In future all insertions should be written to the file
    but only ones in the filtered threshold should be counted in the dictionary.

    """
    # List of tuples of genome features
    genome = fttLookup(organism, experiment)
    # list of tuples of each mapped insertion to be immediately written per insertion
    mappedHitList = []
    # Dictionary with running total of cpm per gene; keys are genes, values are aggregate cpm
    geneDict = {}
    with open('{0}/{1}_bowtie_mapped.txt'.format(experiment, sample), 'r', newline='') as csvfileR:
        sitesReader = csv.reader(csvfileR, delimiter='\t')
        next(sitesReader, None) #skip the headers
        for line in sitesReader:
            contig, nucleotide, Lcounts, Rcounts, totalCounts, cpm = line[0:6]
            nucleotide = int(nucleotide)
            cpm = float(cpm)
            # Used to save previous feature for intergenic calling
            prevFeature = ''
            for gene in genome:
                locus, start, end, strand, length, pid, gene, locus_tag, code, cog, product = \
                    gene[0], int(gene[1]), int(gene[2]), gene[3], gene[4], \
                        gene[5], gene[6], gene[7], gene[8], gene[9], gene[10]
                # number of genes/features hit.
                # 0 = intergenic. 1 = a gene. 2 = overlapping genes.
                # if 2, will print 2 entries in the output.
                genesHit = 0
                # contig from insertion; locus from lookup table
                if contig == locus:
                    # TODO: Add threeprimeness filtering math in here.
                    if nucleotide >= start:
                        if nucleotide <= end:
                            # 0.0 = 5'end ; 1.0 = 3'end
                            # TODO: Should featureEnd have +1 added?
                            if strand == '+':
                                threePrimeness = (nucleotide-start)/(end-start)
                            if strand == '-':
                                threePrimeness = (end-nucleotide)/(end-start)
                            mappedHit = (contig, nucleotide, Lcounts, Rcounts,
                                totalCounts, cpm, threePrimeness, locus_tag)
                            mappedHitList.append(mappedHit)
                            # Add to the total for that gene --
                            # Single-element list (rather than interger) so
                            # that it is subscriptable to add cpm counts
                            # TODO: ThreePrimess math to determine whether
                            # to add to the geneDict
                            geneDict.setdefault(locus_tag, [0])[0] += cpm
                            # Counts for overlap. Not being used currently!
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
                            #        print(genome[locus_tag][i+1])  # the feature after (feature)
                            #    except:
                            #        print(genome[locus_tag][0])
                prevFeature = locus_tag
        genesHit = 0
    # Write individual insertions to *_bowtie_mapped_genes.txt
    with open('{0}/{1}_bowtie_mapped_genes.txt'.format(experiment, sample), 'w', newline='') as csvfileW:
        headers = ('contig', 'nucleotide', 'Lcounts', 'Rcounts', 'totalCounts', 'cpm', 'threePrimeness', 'locus_tag')
        mappedGeneWriter = csv.writer(csvfileW, delimiter='\t')
        mappedGeneWriter.writerow(headers)
        for hit in mappedHitList:
            mappedGeneWriter.writerow(hit)
    # Return aggregated insertions by gene
    return geneDict

def buildGeneTable(organism, sample_list, gene_mappings, experiment=''):
    """
    For each entry in a feature table (.ftt) list the summary of hits
    for each sample in the experiment

    """

    #TODO: Bring back in the header row in the future. Use it here; ignore it for previous steps
    gene_table = fttLookup(organism, experiment)

    # Header will be extended in the future to
    # list the experiment and barcode of each sample of interest
    header = ['Contig', 'Start', 'End', 'Strand', 'Length', 'PID',
        'Gene', 'Synonym', 'Code', 'COG', 'Product']

    # Add header row to ftt file
    gene_table.insert(0, header)

    # Get individual mapped hits
    # mappedHitList = mapToGene(organism, experiment)

    # current column, sample being matched
    currentColumn = len(gene_table[0])-1
    currentSample = ''

    for sampleName in sample_list:
        # Add the new sample name as a new column the table
        currentColumn += 1
        gene_table[0].append(sampleName)
        # Fill the rest of the new column with 0 as the count for each gene
        for row in gene_table[1:]:
            row.append(0)
        # add the sample's results to the building gene_table
        mapped_genes = gene_mappings[sampleName]
        ### mapped_genes = gene_mappings.get(sampleName)
        #TODO: simplify this:
        # for each row in the table, *try* from mapped_genes. Add if found.
        for gene in mapped_genes:
            for i, f in enumerate(gene_table):
                hitLocusTag = gene
                fttLocusTag = f[7]
                # matches based on locusTag.
                # In future should I instead create an index field in the .ftt?
                if hitLocusTag == fttLocusTag:
                    gene_table[i][currentColumn] += mapped_genes[gene][0]
        with open('{0}/summary_gene_table.txt'.format(experiment), 'w') as fo:
            for line in gene_table:
                for entry in line:
                    fo.write('{0}\t'.format(entry))
                fo.write('\n')

def fttLookup(organism, experiment=''):
    """
    Import the ftt file and process as a dictionary of lookup values
    indexed on Synonym (i.e., Locus Tag)
    {'VF_0001': {'locus': 'CP000020', 'start': ...},
        'VF_0002': {'locus': 'CP000020', 'start': ...}}

    """
    # TODO: Error checking when generating the ftt file that locus tags are \
    # unique and complete.
    fttList = []
    with open('{0}/genome_lookup/{1}.ftt'.format(experiment, organism), newline='') as csvfile:
        fttreader = csv.reader(csvfile, delimiter='\t')
        for line in fttreader:
            # ignore header row
            if line[0] != ('Locus'):
                # Locus, Location_Start, Location_End, Strand, Length, PID,
                # Gene, Synonym, Code, COG, Product
                featureData = [line[0], line[1], line[2], line[3], line[4], \
                    line[5], line[6], line[7], line[8], line[9], line[10]]
                fttList.append(featureData)
    return fttList



# ===== Start here ===== #

def main():
    #bowtieOutput = sys.argv[1]
    experiment, organism = 'example01', 'genome'
    #mapSites(bowtieOutput)
    fttLookup(organism, experiment)
    #sample_file = sys.argv[2]
    #experiment = sys.argv[3]
    #processList = samplesToProcess(sample_file, experiment)
    #for s in processList:
    #    print(s)

if __name__ == '__main__':
    main()
