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
        features = False # in the FEATURES section of the GenBank file
        for i, line in enumerate(fi):

            # Don't parse blank lines
            if line.strip():
                parts = line.split()

                # Locus (replicon) as header
                if(parts[0] == 'LOCUS'):
                    locus = parts[1]
                    print('LOCUS\t{}'.format(locus))
                    header = ('Location', 'Strand', 'Length', 'PID',
                        'Gene', 'Synonym', 'Code', 'COG', 'Product')
                    print('\t'.join(header))

                if(parts[0] == 'ORIGIN'):
                    features = False  # Not in FEATURES any more

                if features:
                    if(parts[0] == 'CDS'):

                        # Minus strand if the location begin with 'complement'/'c'
                        if parts[1][0] == 'c':
                            strand = '-'
                            location = parts[1][parts[1].index('(')+1:-1]
                            #location = re.split('..',parts[1])
                        else:
                            strand = '+'
                            location = parts[1]



                        output = (strand, location)
                        print('\t'.join(output))




                if(parts[0] == 'FEATURES'):
                    features = True


                    if i > 5000:
                        break


# ===== Start here ===== #

def main():
    inputfile = 'JNFR01.1.gbff'
    gbk2ftt(inputfile)

if __name__ == "__main__":
    main()
