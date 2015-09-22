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
    parser.add_argument('-g', '--genome',
        help='genome in GenBank format (one concatenated file for multiple contigs/chromosomes)',
        required=True)
    return parser.parse_args(args)

# Change to the specified directory
# http://stackoverflow.com/a/13197763
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)
    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)
    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

# ===== Start here ===== #

def main():
    args = parseArgs(sys.argv[1:])
    gbkfile = args.genome
    experiment = args.experiment
    reads = args.input
    samples = args.samples

    # originally this was to name the organism, e.g. VF for Vfischeri
    # it really only was used to name the genome files so it doesn't
    # seem worth to use it as a variable. Will replace it simply
    # with the word 'genome'.
    # For now, files generated will be called 'genome.fna' etc
    # In the future this can be streamlined.
    organism = 'genome'

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
    readsAssigned = '{0}_assigned.fastq'.format(experiment) # already exsists
    bowtieOutput = '{0}_bowtieOutput.txt'.format(experiment) # will get created at next step

    # Change directory, build bowtie indexes, call bowtie, change directory back
    with cd(tempDir):
        bowtieBuild(organism)
        bowtieMap(organism, readsAssigned, bowtieOutput)

    # Summarize the bowtie results
    # Need more consistency in how directory locations are handled
    insertionNucleotides(tempDir + bowtieOutput, experiment)
    insertionNucleotidesCount(experiment)
    filterSortCounts(experiment)
    #insertionCounts(experiment) ## FOR TESTING ONLY - DOES NOT WRITE DATA
    #normalizeCpm(experiment) ## FOR TESTING ONLY UNLESS WE WANT IT TO WRITE DATA
    #mapToGene(organism, experiment) ## FOR TESTING ONLY UNLESS WE WANT IT TO WRITE DATA
    mapToGeneSummary(1.0, organism, experiment)  # First parameter is threePrimeness cutoff

if __name__ == '__main__':
    main()
