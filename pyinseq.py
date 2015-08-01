#!/usr/bin/env python
"""
Main script for running the pyinseq package

"""

from assign import *
from gbkconvert import *
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


# ===== Start here ===== #

def main():
    #args = parse_args(sys.argv[1:])
    #print(args.input)
    #print(args.samples)
    #assign_and_trim(args.input, args.samples)
    inputfile = sys.argv[1]
    organism = sys.argv[2]
    experiment = sys.argv[3]

    # Create the directory struture based on the experiment name
    create_directories(experiment)

    # Prepare genome files from the GenBank input
    gbk2fna(inputfile, organism, experiment)
    gbk2ftt(inputfile, organism, experiment)

    # Prepare the bowtie indexes
    


if __name__ == '__main__':
    main()
