#!/usr/bin/env python
"""

GenBank conversion utilities for the pyinseq pipeline.

# gbk2fasta()
Convert GenBank sequence to fasta sequence
Multilocus GenBank converts to one multifasta GenBank file
Locus headers are the fasta headers

# gbk2ftt()
Convert GenBank to feature table
Format of .ptt and .rnt files, including the following features:
    CDS
    rRNA
    tRNA
    misc_RNA
Multilocus GenBank converts to multiple .ftt files
Locus headers are the file names (locus.ftt)
Unlike .ptt files that show the number of amino acids as 'length'

(BETTER TO HAVE ONE FILE -- RIGHT?)

"""

import re

def gbk2fasta(infile):
    with open(infile, 'r') as fi:
        dna_seq = False # in the DNA sequence of the file
        for i, line in enumerate(fi):

            # Don't parse blank lines
            if line.strip():
                parts = line.split()

                # Locus (replicon) as header
                if(parts[0] == 'LOCUS'):
                    locus = parts[1]
                    print('>{}'.format(locus))

                # DNA Sequence
                if(parts[0] == '//'):
                    dna_seq = False
                if dna_seq:
                    sequence = ''.join(n for n in line.strip() if n.isalpha())
                    print(sequence)
                if(parts[0] == 'ORIGIN'):
                    dna_seq = True

def gbk2ftt(infile):
    with open(infile, 'r') as fi:

        # Initialize variables
        features = False # in the FEATURES section of the GenBank file
        new_feature = False # feature to be written
        location = '0..0'
        strand = '+'
        length = 0
        protein_id = ''
        gene = '-'
        locus_tag = ''
        cog = '-'
        product = ''
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
                    print('LOCUS\t{}'.format(locus))
                    header = ('Location', 'Strand', 'Length', 'PID',
                        'Gene', 'Synonym', 'Code', 'COG', 'Product')
                    print('\t'.join(header))

                if(parts[0] == 'ORIGIN'):
                    features = False  # Not in FEATURES any more


                if features:

                    # Only features CDS, etc. here will flip the
                    # new_feature to True and will be written in the ouput

                    if(parts[0] == 'CDS'):
                        new_feature = True
                        # Remove any < or > characters
                        # Minus strand if the location begin with 'complement'/'c'
                        if parts[1][0] == 'c':
                            strand = '-'
                            location = parts[1][parts[1].index('(')+1:-1].strip('<>')
                        else:
                            strand = '+'
                            location = parts[1].strip('<>')

                        # Recall stripped out < and > so length may not be
                        # divisible by 3 for CDS if gene is at end of contig.
                        first = location[0:location.index('..')]
                        last = location[location.index('..')+2:]
                        length = int(last )- int(first) + 1

                    if '/protein_id=' in parts[0]:
                        protein_id = parts[0][13:-1]

                    if '/gene=' in parts[0]:
                        gene = parts[0][7:-1]

                    if '/locus_tag=' in parts[0]:
                        locus_tag = parts[0][12:-1]

                    if product_append:
                        # /product extends > 2 lines
                        if not line.strip()[-1] == '\"':
                            product = product + ' ' + line.strip()
                            product_append = True
                        # This is last line of the /product field
                        else:
                            product = product + ' ' + line.strip()[:-1]
                            product_append = False

                    if '/product=' in parts[0]:
                        product = line.strip()[10:-1]
                        if not line.strip()[-1] == '\"':
                            product = line.strip()[10:]
                            product_append = True

                    # Print line if all data are present (gene, COG not required)
                    # Reset flags/defaults
                    if new_feature:
                        if locus_tag:
                            if protein_id:
                                if not product_append:
                                    output = (location, strand, str(length),
                                        protein_id, gene, locus_tag, cog, product)
                                    print('\t'.join(output))
                                    new_feature = False
                                    protein_id = ''
                                    gene = '-'
                                    locus_tag = ''
                                    cog = '-'
                                    product = ''

                if(parts[0] == 'FEATURES'):
                    features = True


                    if i > 10000:
                        break


# ===== Start here ===== #

def main():
    inputfile = 'JNFR01.1.gbff'
    gbk2ftt(inputfile)

if __name__ == "__main__":
    main()
