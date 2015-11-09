#!/usr/bin/env python
"""
Main script for running the pyinseq package

"""

#from assign import *
from demultiplex import *
from gbkconvert import *
from mapReads import *
from processMapping import *
from utils import *
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
    parser.add_argument('-d', '--disruption',
        help='fraction of gene disrupted (0.0 - 1.0)',
        default=1.0)
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
    experiment = convert_to_filename(args.experiment)
    reads = args.input
    samples = args.samples
    disruption = float(args.disruption)

    # originally this was to name the organism, e.g. VF for Vfischeri
    # it really only was used to name the genome files so it doesn't
    # seem worth to use it as a variable. Will replace it simply
    # with the word 'genome'.
    # For now, files generated will be called 'genome.fna' etc
    # In the future this can be streamlined.
    organism = 'genome'

    # pyinseqDirectory = os.getcwd()
    genomeDir = '{experiment}/genome_lookup/'.format(experiment=experiment)
    # Create the directory struture based on the experiment name
    createExperimentDirectories(experiment)

    #demultiplex based on barcodes defined in the sample file
    demultiplex_fastq(reads, samples, experiment)

    # Prepare genome files from the GenBank input
    gbk2fna(gbkfile, organism, genomeDir)
    gbk2ftt(gbkfile, organism, genomeDir)

    # Change directory, build bowtie indexes, change directory back
    with cd(genomeDir):
        bowtieBuild(organism)

    # List of file paths of the demultiplexed files
    demultiplexedSample_list = demultiplexedSamplesToProcess(samples, experiment)

    # PROCESS SAMPLE
    # Rewrite as new .fastq file with only chromosome sequence
    # Trim barcode, trim transposon. Trim corresponding quality.
    # Ignore if not a good transposon junction.
    for samplePath in demultiplexedSample_list:
        #TODO: Fix so that the sample name does not need to be calculated
        s1 = samplePath.split('.')[0].rfind('/')
        sampleName = samplePath.split('.')[0][s1+1:]
        trimmedSamplePath = '{experiment}/{sampleName}_trimmed.fastq'.format(
            experiment=experiment,
            sampleName=sampleName)
        trim_fastq(samplePath, trimmedSamplePath, sampleName)
        # Change directory, map to bowtie, change directory back
        trimmedSampleFile = '{0}_trimmed.fastq'.format(sampleName)
        bowtieOutputFile = '{0}_output_bowtie.txt'.format(sampleName)
        with cd(genomeDir):
            # Paths are relative to the genome_lookup directory
            # from where bowtie is called
            bowtie_in = '../{0}'.format(trimmedSampleFile)
            bowtie_out = '../{0}'.format(bowtieOutputFile)
            bowtieMap(organism, bowtie_in, bowtie_out)
        mapSites('{0}/{1}'.format(experiment, bowtieOutputFile))



    # Assign and trim barcodes
    # Now currently filtering by default (16-17 bp, TA at end)
    # In future could add as command line option
    #assignAndTrim(reads, samples, experiment, tempDir)


    # Summarize the bowtie results
    # Need more consistency in how directory locations are handled
#    insertionNucleotides(tempDir + bowtieOutput, experiment)
#    insertionNucleotidesCount(experiment)
#    filterSortCounts(experiment)
    #insertionCounts(experiment) ## FOR TESTING ONLY - DOES NOT WRITE DATA
    #normalizeCpm(experiment) ## FOR TESTING ONLY UNLESS WE WANT IT TO WRITE DATA
    #mapToGene(organism, experiment) ## FOR TESTING ONLY UNLESS WE WANT IT TO WRITE DATA
#    mapToGeneSummary(disruption, organism, experiment)

if __name__ == '__main__':
    main()
