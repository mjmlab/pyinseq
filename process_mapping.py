#!/usr/bin/env python
"""
Counts the bowtie hits at each position in each sample

Future - filter on the 16/17 bp positions before this step; maybe even before bowtie mapping

"""

import sys
import re
from collections import Counter

def multifasta2dict(fna):
    """
    Returns a dictionary of fasta header and sequences

    Input file in multifasta format
    >header1
    agctag
    >header2
    gattta

    Function returns the fasta data as a dictionary
    {'header1': 'agctag', 'header2': 'gattta'}

    """
    sequence_dict = {}
    with open(fna, 'r') as fi:
        header = ''
        first_line = False # next line is first line of sequence in contig
        for line in fi:
            line = line.rstrip()
            if line.startswith('>'): # header line; therefore new contig
                header = line[1:] # new header
                first_line = True
            else: # sequence line
                if first_line:
                    sequence_dict[header] = line # create dictionary entry
                    first_line = False
                # would this work in one step?
                else:
                    appended_sequence = sequence_dict[header] + line
                    sequence_dict[header] = appended_sequence
        print(sequence_dict)
        return sequence_dict

def TA_sites(fna):
    """ List TA dinucleotides in a fasta nucleotide file

    Returns a list of nucleotide positions for
    5'-TA for each contig in the fna file
    e.g., CP000020 = (20, 22, 24, ... 1401104)
    #INSERT THE REAL DATA HERE! HOW DOES IT WORK FOR MULTIFASTA?
    Assumes circular genome
    i.e. will check if the last nucleotide is a T and the first is an A
    """
    with open(fna, 'r') as fi:
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


        contig_length = len(fna_string)

        print(contig_name)
        print(len(fna_string))

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

def normalize_cpm(bowtie_output):
    """ insert """
    pass

def map_to_gene(normalized_output):
    """ insert """
    pass


# ===== Start here ===== #

def main():
    fna = sys.argv[1]
    multifasta2dict(fna)

if __name__ == "__main__":
    main()
