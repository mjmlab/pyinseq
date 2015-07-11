#!/usr/bin/env python
"""
Counts the bowtie hits at each position in each sample

Future - filter on the 16/17 bp positions before this step; maybe even before bowtie mapping

"""

import sys
import re
from collections import Counter

def count_bowtie(bowtie_output):
    """ List counts of reads to the Left and Right of each TA site

    future:
    default show_all=False
    show_all=True will show all of the sites even those will count=0

    THIS MODULE STILL NEEDS TO BE WRITTEN!!!!!!!!

    """
    with open(bowtie_output, 'r') as fi:
        outfile = bowtie_output + '_counted'
        with open(outfile, 'w') as fo:
            count = {}
            for line in fi:
                sample_assignment = line.split('\t')[0].rfind('//')
                read_data_to_count = ('{0}\t{1}\t{2}\t{3}'.format(
                    line.split('\t')[0][sample_assignment+2:],
                    line.split('\t')[1],
                    line.split('\t')[2],
                    line.split('\t')[3]))
                count[read_data_to_count] = count.get(read_data_to_count, 0) + 1
            for c in count:
                print('{0}\t{1}'.format(c, count[c]))

def TA_sites(fna):
    """ List TA dinucleotides in a fasta nucleotide file

    Returns a list of TA sites for each contig in the fna file
    e.g., CP000020 = (20, 22, 24, ... 1401104)
    #INSERT THE REAL DATA HERE! HOW DOES IT WORK FOR MULTIFASTA?
    Assumes circular genome
    i.e. will check if the last nucleotide is a T and the first is an A
    """
    with open(fna, 'r') as fi:
        outfile = fna + '_ta_sites.txt'
        with open(outfile, 'w') as fo:
            count = 0
            # need to deal with last base in contig
            fna_string = ''
            contig_name = ''
            for line in fi:
                if line.startswith('>'):
                    if len(fna_string) > 0:
                        # do something with the old DNA sequence before
                        # starting with the new contig
                        pass
                    contig_name = line.rstrip('\n')[1:]
                    fna_string = ''  # New contig...clear the DNA sequence
                else:
                    fna_string = fna_string + line.strip('\n')

            # START WORKING HERE
            ta_list = []
            i = 0   # nucleotide index
            # Checks the linear molecule (i.e. all except the last base)
            while i < len(fna_string)-1:
                if fna_string[i:i+2] == 'ta':
                    ta_list.append(i+1)
                i += 1
            # Checks last base circular around to first base
            if i == len(fna_string)-1:
                if (fna_string[-1] + fna_string[0]) == 'ta':
                    ta_list.append(i+1)
            print(ta_list)

            # NEED TO RENAME IT AS THE CONTIG NAME.
            # STILL NEED TO FIGURE OUT MULTIFASTA


            contig_length = len(fna_string)

            print(contig_name)
            print(len(fna_string))

            # How best to return the list of lists (for a multifasta)?


            """for n in line.rstrip():
                print(len(n))
                fna_list.append(n)
                if(len(fna_list) == 2):
                    # special case of last base matching
                    # need robust way to count number of bases
                    # and print its position here
                    pass
                else:
                    if fna_list[-2:] == ['t','a']:
                        print(len(fna_list)-2)
                if len(fna_list) > 100:
                    exit(0)"""



def map_to_gene(bowtie_output):
    """ insert """
    pass

def map_to_gene(bowtie_output):
    """ insert """
    pass


# ===== Start here ===== #

def main():
    fna = sys.argv[1]
    TA_sites(fna)

if __name__ == "__main__":
    main()
