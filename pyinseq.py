#!/usr/bin/env python
"""
Main script for running the pyinseq package

"""

from assign import *
from gbkconvert import *
from mapReads import *
from processMapping import *
import sys
import os
import argparse

def parseArgs(args):
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
    return parser.parseArgs(args)

# Change to the genome directory and run bowtieBuild() and bowtieMap()
class cdCallBowtie:
    def __init__(self, tempDir, organism, reads, bowtieOutput):
        self.savedPath = os.getcwd()
        os.chdir(tempDir)
        bowtieBuild(organism)
        bowtieMap(organism, reads, bowtieOutput)
    def __del__(self):
        os.chdir(self.savedPath)

# ===== Start here ===== #

def main():
    #args = parseArgs(sys.argv[1:])
    #print(args.input)
    #print(args.samples)
    gbkfile = sys.argv[1]
    organism = sys.argv[2]
    experiment = sys.argv[3]
    reads = sys.argv[4]
    samples = sys.argv[5]

    #pyinseqDirectory = os.getcwd()
    tempDir = '{0}/temp/'.format(experiment)

    # Create the directory struture based on the experiment name
    createDirectories(experiment)

    # Prepare genome files from the GenBank input
    gbk2fna(gbkfile, organism, tempDir)
    gbk2ftt(gbkfile, organism, tempDir)

    # Assign and trim barcodes
    # Now currently filtering by default (16-17 bp, TA at end)
    # In future could add as command line option
    assignAndTrim(reads, samples, experiment, tempDir)

    # Prepare the bowtie indexes
    # Map the reads using bowtie_map
    readsAssigned = '{0}_Assigned.fastq'.format(experiment) # already exsists
    bowtieOutput = '{0}_BowtieOutput.txt'.format(experiment) # will get created at next step
    x = cdCallBowtie(tempDir, organism, readsAssigned, bowtieOutput)



if __name__ == '__main__':
    main()
