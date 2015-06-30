#!/usr/bin/env python
"""
Counts the bowtie hits at each position in each sample

"""

import sys
import re
from collections import Counter

def count_bowtie(bowtie_output):

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

# ===== Start here ===== #

def main():
    bowtie_output = sys.argv[1]
    count_bowtie(bowtie_output)

if __name__ == "__main__":
    main()
