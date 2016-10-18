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

import csv
import os
import re
import sys


def gbk2fna(infile, organism, outputdirectory=''):
    """Convert genbank format to fna format."""
    with open(infile, 'r') as fi:
        outfile = '{0}{1}.fna'.format(outputdirectory, organism)
        print('  Nucleotide file output file: {}'.format(outfile))
        if not os.path.exists('{0}'.format(outputdirectory)):
            print('Error: {0} directory was not created.'.format(outputdirectory))
            exit(1)
        with open(outfile, 'w') as fo:
            dna_seq = False  # in the DNA sequence of the file
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
    """Convert genbank format to ptt-like ftt format."""
    with open(infile, 'r') as fi:
        outfile = '{0}{1}.ftt'.format(outputdirectory, organism)
        print('  Feature table output file: {}'.format(outfile))
        with open(outfile, 'w') as fo:
            writer = csv.writer(fo, delimiter='\t', dialect='excel')
            header = ('Locus', 'Location_Start', 'Location_End', 'Strand', 'Length', 'PID',
                      'Gene', 'Synonym', 'Code', 'COG', 'Product')
            writer.writerow(header)

            # Initialize variables
            features = False  # in the FEATURES section of the GenBank file
            new_feature = False  # collecting data for a new feature
            parse_types = ['CDS', 'tRNA', 'rRNA', 'misc_RNA']
            strand = '+'
            length = 0
            protein_id = '-'
            gene = '-'
            locus_tag = '-'
            code = '-'
            cog = '-'
            product = '-'
            product_append = False  # append the current line to product

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
                                # if locus_tag:
                                if not product_append:
                                    output = (locus, first, last, strand,
                                              str(length), protein_id, gene,
                                              locus_tag, code, cog, product)
                                    writer.writerow(output)
                                    new_feature = False

                        if line[5:21].rstrip() in parse_types:
                            new_feature = True  # Feature that should be written
                            protein_id = '-'
                            gene = '-'
                            locus_tag = '-'
                            code = '-'
                            cog = '-'
                            product = '-'

                            # NOTES ABOUT FEATURES
                            # 1. At ends of contigs greater than/less than signs
                            #    (> / <) are removed.
                            # 2. Complicated features use only the outer bounds
                            #    join(481257..481331,481333..482355) uses 481257..482355
                            location = re.search(r'(\d+)\.+.*\.(\d+)', parts[1])
                            first = location.group(1)
                            last = location.group(2)
                            try:
                                location_raw = '{0}..{1}'.format(first, last)
                                length = int(last) - int(first) + 1
                            except AttributeError:
                                errorComplexFeature = \
                                'PyINSeq Error: Complex feature coordinates at or near {0} ' \
                                'in GenBank file. Additional attention required.'.format(locus_tag)
                                print(errorComplexFeature)
                                exit(0)
                            strand = '-' if parts[1].startswith('complement') else '+'

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
                            # TODO(Error handling if not exactly 2 parentheses)
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


def main():
    """Start here."""
    inputfile = sys.argv[1]
    organism = sys.argv[2]
    # gbk2fna(inputfile, organism)
    gbk2ftt(inputfile, organism)

if __name__ == '__main__':
    main()
