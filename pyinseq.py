#!/usr/bin/env python
"""
Main script for running the pyinseq package

"""

from assign import *
from gbkconvert import *
from map_reads import *
import sys
import os
import argparse

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',
        help='input Illumina reads file',
        required=True)
    parser.add_argument('-s', '--samples',
        help='sample list with barcodes',
        required=True)
    parser.add_argument('-e', '--experiment',
        help='experiment name (no spaces or special characters)',
        required=True)
    return parser.parse_args(args)

# Change to the genome directory and run bowtie_build() and bowtie_map()
class cd_bowtie:
    def __init__(self, temp_dir, organism, reads, bowtie_output):
        self.savedPath = os.getcwd()
        os.chdir(temp_dir)
        bowtie_build(organism)
        bowtie_map(organism, reads, bowtie_output)
    def __del__(self):
        os.chdir(self.savedPath)

# ===== Start here ===== #

def main():
    #args = parse_args(sys.argv[1:])
    #print(args.input)
    #print(args.samples)
    #assign_and_trim(args.input, args.samples)
    gbkfile = sys.argv[1]
    organism = sys.argv[2]
    experiment = sys.argv[3]
    reads = sys.argv[4]
    samples = sys.argv[5]

    #pyinseq_directory = os.getcwd()
    temp_dir = '{0}/temp/'.format(experiment)


    # Create the directory struture based on the experiment name
    create_directories(experiment)

    # Prepare genome files from the GenBank input
    gbk2fna(gbkfile, organism, temp_dir)
    gbk2ftt(gbkfile, organism, temp_dir)

    # Assign and trim barcodes
    assign_and_trim(reads, samples, experiment, temp_dir)



    # Prepare the bowtie indexes
    # Map the reads using bowtie_map
    #x = cd_bowtie(temp_dir, organism, reads, '../bowtie_output.txt')


if __name__ == '__main__':
    main()
