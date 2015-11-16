#!/usr/bin/env python
"""

GenBank conversion utilities for the pyinseq pipeline.

# gbk2fna()
Convert GenBank sequence to fasta sequence
Multilocus GenBank converts to one multifasta GenBank file
Locus headers are the fasta headers
Maintains original newlines (typically leaving up to 60 nucleotides per line)
File is written to temp/ directory for the EXPERIMENT:
    EXPERIMENT/temp/genome.fna

# gbk2ftt()
Convert GenBank to feature table
Format similar to .ptt and .rnt files except:
- full tabular (locus as a field)
- start and end positions as separate fields
- includes the following features:
    CDS
    rRNA
    tRNA
    misc_RNA
- Multilocus GenBank converts to multi-.ftt file
'Unlike .ptt files that show the number of amino acids as 'length'

"""

import sys
import os
import re
import csv

def gbk2fna(infile, organism, outputdirectory=''):
    with open(infile, 'r') as fi:
        outfile = '{0}{1}.fna'.format(outputdirectory, organism)
        if not os.path.exists('{0}'.format(outputdirectory)):
            print('Error: {0} directory was not created.'.format(outputdirectory))
            exit(1)
        with open(outfile, 'w') as fo:
            dna_seq = False # in the DNA sequence of the file
            for i, line in enumerate(fi):

                # Don't parse blank lines
                if line.strip():
                    parts = line.split()

                    # Locus (replicon) as header
                    if(parts[0] == 'LOCUS'):
                        locus = parts[1]
                        fo.write('>{}\n'.format(locus))

                    # DNA Sequence
                    if(parts[0] == '//'):
                        dna_seq = False
                    if dna_seq:
                        sequence = ''.join(n for n in line.strip() if n.isalpha())
                        fo.write('{0}\n'.format(sequence))
                    if(parts[0] == 'ORIGIN'):
                        dna_seq = True

def gbk2ftt(infile, organism, outputdirectory=''):
    with open(infile, 'r') as fi:
        outfile = '{0}{1}.ftt'.format(outputdirectory, organism)
        with open(outfile, 'w') as fo:
            writer = csv.writer(fo, delimiter='\t', dialect='excel')
            header = ('Locus', 'Location_Start', 'Location_End', 'Strand', 'Length', 'PID',
                'Gene', 'Synonym', 'Code', 'COG', 'Product')
            writer.writerow(header)

            # Initialize variables
            features = False # in the FEATURES section of the GenBank file
            new_feature = False # collecting data for a new feature
            new_feature_type = ''
            parse_types = ['CDS', 'tRNA', 'rRNA', 'misc_RNA']
            location = '0..0'
            strand = '+'
            length = 0
            protein_id = '-'
            gene = '-'
            locus_tag = '-'
            code = '-'
            cog = '-'
            product = '-'
            product_append = False # append the current line to product

            for i, line in enumerate(fi):

                # Don't parse blank lines
                if line.strip():
                    parts = line.split()

                    # PRINT HEADER FOR THE LOCUS
                    # 2 LINES:
                    # LOCUS <tab> locus name
                    # Location <tab> Strand etc...
                    if(parts[0] == 'LOCUS'):
                        locus = parts[1]

                    if features:
                        # Print line before go on to next feature
                        # (gene, COG, protein id not required)
                        # Reset flags/defaults
                        if line[5:21].rstrip():
                            if new_feature:
                                #if locus_tag:
                                if not product_append:
                                    output = (locus, first, last,
                                        strand, str(length), protein_id, gene,
                                        locus_tag, code, cog, product)
                                    writer.writerow(output)
                                    new_feature = False
                                    new_feature_type = ''
                                    protein_id = '-'
                                    gene = '-'
                                    locus_tag = '-'
                                    code = '-'
                                    cog = '-'
                                    product = '-'

                        if(line[5:21].rstrip() in parse_types):
                            new_feature = True # Feature that should be written
                            new_feature_type = line[5:21].rstrip()


                            # SIMPLE FEATURES - TWO COORDINATES, FORWARD OR COMPLEMENT
                            # Minus strand if the location begin with 'complement'/'c'
                            if parts[1][0] == 'c':
                                strand = '-'
                                location_raw = parts[1][parts[1].index('(')+1:-1]
                            else:
                                strand = '+'
                                location_raw = parts[1]


                            # MILDLY COMPLICATED FEATURES
                            # e.g., join(481257..481331,481333..482355)
                            # it just reports outer bounds: 481257..482355
                            # assumes not too complicated! - feature on same strand, on same contig

                            # Regex to pull out start location = (?<=\()(\d+)(?=\.\.)
                                # Matches digits in: (digits..
                            # Regex to pull out end location = (?<=\.\.)(\d+)(?=\))
                                # Matches digits in: ..digits)

                            # Addressses this case: 'join(481257..481331,481333..482355)'
                            # TODO: join(complement(1..5,7..10))
                            if parts[1].startswith('join'):
                                start_match = re.search(r'(?<=\()(\d+)(?=\.\.)', parts[1])
                                end_match = re.search(r'(?<=\.\.)(\d+)(?=\))', parts[1])
                                try:
                                    location_raw = '{0}..{1}'.format(start_match.group(), end_match.group())
                                except AttributeError:
                                    errorComplexFeature = \
                                    'PyINSeq Error: Complex feature coordinates at or near {0} ' \
                                    'in GenBank file. Additional attention required.'.format(locus_tag)
                                    print(errorComplexFeature)
                                    exit(0)

                            # Process location info:
                            # Strip out < and > and note that may not be
                            # ... divisible by 3 for CDS if gene is at end of contig
                            first = location_raw[:location_raw.find('..')].strip('<>')
                            last = location_raw[location_raw.rfind('..')+2:].strip('<>')
                            location = '{0}..{1}'.format(first, last)
                            length = int(last) - int(first) + 1

                        if '/protein_id=' in parts[0]:
                            protein_id = parts[0][13:-1]

                        if '/gene=' in parts[0]:
                            gene = parts[0][7:-1]

                        if '/locus_tag=' in parts[0]:
                            locus_tag = parts[0][12:-1]

                        # Multi-line product description
                        if product_append:
                            product = product + ' ' + line.strip()
                            if product.count('\"') != 2:
                                product_append = True
                            # TODO: Error handling if not exactly 2 parentheses
                            if product.count('\"') == 2:
                                product = product.strip('\"')
                                product_append = False

                        if '/product=' in parts[0]:
                            product = line.strip()[9:]
                            if product.count('\"') != 2:
                                product_append = True
                            if product.count('\"') == 2:
                                product = product.strip('\"')

                    if(parts[0] == 'ORIGIN'):
                        features = False  # Not in FEATURES any more
                    if(parts[0] == 'FEATURES'):
                        features = True

# ===== Start here ===== #

def main():
    inputfile = sys.argv[1]
    organism = sys.argv[2]
    #gbk2fna(inputfile, organism)
    gbk2ftt(inputfile, organism)

if __name__ == '__main__':
    main()
