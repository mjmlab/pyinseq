#!/usr/bin/env python
"""

GenBank conversion utilities for the pyinseq pipeline.

# gbk2fna()
Convert GenBank sequence to fasta sequence
Multilocus GenBank converts to one multifasta GenBank file
Locus headers are the fasta headers
Maintains original newlines (typically leaving up to 60 nucleotides per line)

# gbk2ftt()
Convert GenBank to feature table
Format similar to .ptt and .rnt files, including the following features:
    CDS
    rRNA
    tRNA
    misc_RNA
Multilocus GenBank converts to multi-.ftt file
Locus headers are the file names (locus.ftt)
Unlike .ptt files that show the number of amino acids as 'length'

Future:
Update with argument passing
Add argument passing
UPDATE INFO ABOVE. INFO ON FILE NAMES, ETC.
Extract COG from /note field in Refseq GBK files
Allow user to generate ptt files by filtering to CDS and dividing length by 3
Or CDS + misc_RNA but exclude rRNA/tRNA
Consider parsing on fixed width rather than on space-delimited for LOCUS, etc

"""

import sys

def gbk2fna(infile, organism):
    with open(infile, 'r') as fi:
        outfile = organism + '.fna'
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

def gbk2ftt(infile, organism):
    with open(infile, 'r') as fi:
        outfile = organism + '.ftt'
        with open(outfile, 'w') as fo:

            # Initialize variables
            features = False # in the FEATURES section of the GenBank file
            new_feature = False # collecting data for a new feature
            new_feature_type = ''
            parse_types = ['CDS', 'tRNA', 'rRNA', 'misc_RNA']
            location = '0..0'
            strand = '+'
            length = 0
            protein_id = ''
            gene = '-'
            locus_tag = ''
            code = '-'
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
                        fo.write('LOCUS\t{}\n'.format(locus))
                        header = ('Location', 'Strand', 'Length', 'PID',
                            'Gene', 'Synonym', 'Code', 'COG', 'Product')
                        fo.write('\t'.join(header) + '\n')

                    if features:

                        # Print line before go on to next feature
                        # (gene, COG, protein id not required)
                        # Reset flags/defaults
                        if line[5:21].rstrip():
                            if new_feature:
                                #if locus_tag:
                                if not product_append:
                                    output = (location, strand, str(length),
                                        protein_id, gene, locus_tag,
                                        code, cog, product)
                                    fo.write('\t'.join(output) + '\n')
                                    new_feature = False
                                    new_feature_type = ''
                                    protein_id = ''
                                    gene = '-'
                                    locus_tag = ''
                                    code = '-'
                                    cog = '-'
                                    product = ''

                        if(line[5:21].rstrip() in parse_types):
                            new_feature = True # Feature that should be written
                            new_feature_type = line[5:21].rstrip()

                            # Remove any < or > characters
                            # Minus strand if the location begin with 'complement'/'c'
                            if parts[1][0] == 'c':
                                strand = '-'
                                location_raw = parts[1][parts[1].index('(')+1:-1]
                            else:
                                strand = '+'
                                location_raw = parts[1]

                            # Process location info:
                            # Strip out < and > and note that may not be
                            # ... divisible by 3 for CDS if gene is at end of contig
                            first = location_raw[0:location_raw.index('..')].strip('<>')
                            last = location_raw[location_raw.index('..')+2:].strip('<>')
                            location = '{0}..{1}'.format(first, last)
                            length = int(last)- int(first) + 1

                        if '/protein_id=' in parts[0]:
                            protein_id = parts[0][13:-1]

                        if '/gene=' in parts[0]:
                            gene = parts[0][7:-1]

                        if '/locus_tag=' in parts[0]:
                            locus_tag = parts[0][12:-1]

                        #if product_append:
                            # /product extends > 2 lines
                            #if not line.strip()[-1] == '\"':
                                #product = product + ' ' + line.strip()
                                #product_append = True
                            # This is last line of the /product field
                            #else:
                                #product = product + ' ' + line.strip()[:-1]
                                #product_append = False"""
                        if '/product=' in parts[0]:
                            product = line.strip()[10:-1]
                            #if not line.strip()[-1] == '\"':
                                #product = line.strip()[10:]
                                #product_append = True

                    if(parts[0] == 'ORIGIN'):
                        features = False  # Not in FEATURES any more
                    if(parts[0] == 'FEATURES'):
                        features = True

                        #if i > 5000:
                            #break


# ===== Start here ===== #

def main():
    inputfile = sys.argv[1]
    organism = sys.argv[2]
    #parser = argparse.ArgumentParser()
    #parser.add_argument('-f', '--fasta', action="store_true", default=False)
    #parser.add_argument('-t', '--table', action="store_true", default=False)

    gbk2fna(inputfile, organism)
    gbk2ftt(inputfile, organism)

if __name__ == "__main__":
    main()
