#!/usr/bin/env python
"""

GenBank conversion utilities for the pyinseq pipeline.

# gbk2fasta()
Convert GenBank sequence to fasta sequence
Multilocus GenBank converts to one multifasta GenBank file
Locus headers are the fasta headers

# gbk2ftt()
Convert GenBank to feature table
Format of .ptt and .rnt files, including protein and RNA genes
Multilocus GenBank converts to multiple .ftt files
Locus headers are the file names (locus.ftt)

(BETTER TO HAVE ONE FILE -- RIGHT?)

"""

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




# ===== Start here ===== #

def main():
    inputfile = 'JNFR01.1.gbff'
    gbk2ftt(inputfile)

if __name__ == "__main__":
    main()
