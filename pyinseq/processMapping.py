#!/usr/bin/env python3
'''
Counts the bowtie hits at each position in each sample

'''

import csv
import os

def map_sites(sample, samplesDict, settings):
    '''Map insertions to nucleotide sites.'''
    # Placeholder for dictionary of mapped reads in format:
    # {(contig, position) : [Lcount, Rcount]}
    mapDict = {}
    # overallTotal = denominator for cpm calculation
    overallTotal = 0
    cpm = 0
    bowtie_file = settings.path + sample + '_bowtie.txt'
    sites_file = settings.path + sample + '_sites.txt'
    with open(bowtie_file, 'r') as fi:
        for line in fi:
            bowtiedata = line.rstrip().split('\t')
            # Calculate transposon insertion point = transposonNT
            contig, insertionNT, readLength = str(bowtiedata[2]), int(bowtiedata[3]), len(bowtiedata[4])
            if bowtiedata[1] == '+': # positive strand read
                insertionNt = insertionNT + readLength - 1
                mapDict.setdefault((contig,insertionNt),[0,0])[0] += 1   # Lcount
            else: # negative strand read
                insertionNt = insertionNT + 1
                mapDict.setdefault((contig,insertionNt),[0,0])[1] += 1   # Rcount
            overallTotal += 1
    # write tab-delimited of contig/nucleotide/Lcount/Rcount/TotalCount/cpm
    # use the index totalCounts as the denominator for cpm calculation
    with open(sites_file, 'a') as fo:
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
    return mapDict

def map_genes(organism, sample, disruption, settings):
    '''Maps insertions to genes

       1. Writes a csv file listing the gene for each insertion with the tabs:
       contig, nucleotide, Lcounts, Rcounts, totalCounts, cpm, threePrimeness, locus_tag
       2. Returns a dictionary of aggregate counts per gene (for any gene with at
       least one hit)

       Insertions that map to multiple genes result in multiple lines in the CSV file
       and are counted for both genes in the returned dictionary.

       ThreePrimeness = insertion location in gene (5' end = 0.0, 3' end = 1.0)
       All insertions are written to the file but only ones <= disruption threshold
       are counted in the dictionary.
    '''

    # List of tuples of genome features
    genome = fttLookup(organism, settings.experiment)
    # list of tuples of each mapped insertion to be immediately written per insertion
    mappedHitList = []
    # Dictionary with running total of cpm per gene; keys are genes, values are aggregate cpm
    # only hits in the first part of the gene are added to the count, as defined
    # by the disruption threshold.
    # if disruption = 1.0 then every hit in the gene is included
    geneDict = {}
    sites_file = settings.path + sample + '_sites.txt'
    genes_file = settings.path + sample + '_genes.txt'
    with open(sites_file, 'r', newline='') as csvfileR:
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
                # contig from insertion; locus from lookup table
                if contig == locus:
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
                            # Filter based on location in the gene
                            if threePrimeness <= disruption:
                                # Add to the total for that gene --
                                # Single-element list (rather than interger) so
                                # that it is subscriptable to add cpm counts
                                geneDict.setdefault(locus_tag, [0])[0] += cpm
                prevFeature = locus_tag
    # Write individual insertions to *_genes.txt
    with open(genes_file, 'w', newline='') as csvfileW:
        headers = ('contig', 'nucleotide', 'left_counts', 'right_counts', 'total_counts', 'cpm', 'three_primeness', 'locus_tag')
        mappedGeneWriter = csv.writer(csvfileW, delimiter='\t')
        mappedGeneWriter.writerow(headers)
        for hit in mappedHitList:
            mappedGeneWriter.writerow(hit)
    # Return aggregated insertions by gene (filtered on 5'-3' threshold)
    return geneDict

def build_gene_table(organism, sample_dict, gene_mappings, experiment=''):
    '''
       For each entry in a feature table (.ftt) list the summary of hits
       for each sample in the experiment
    '''

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

    for sampleName in sample_dict:
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
        with open('results/{0}/summary_gene_table.txt'.format(experiment), 'w') as fo:
            writer = csv.writer(fo, delimiter='\t', dialect='excel')
            writer.writerows(gene_table)


def fttLookup(organism, experiment=''):
    '''Import the ftt file and process as a dictionary of lookup values
       indexed on Synonym (i.e., Locus Tag)
       {'VF_0001': {'locus': 'CP000020', 'start': ...},
           'VF_0002': {'locus': 'CP000020', 'start': ...}}
    '''

    # TODO: Error checking when generating the ftt file that locus tags are \
    # unique and complete.
    fttList = []
    with open('results/{0}/genome_lookup/{1}.ftt'.format(experiment, organism), newline='') as csvfile:
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


def main():
    '''Start here.'''
    pass

if __name__ == '__main__':
    main()
