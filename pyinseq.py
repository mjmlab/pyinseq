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

    # samples dictionary
    # samples = OrderedDict([('name1', {'name': 'name1', 'barcode': 'barcode1'}),
    #    ('name2', {'name': 'name2', 'barcode': 'barcode2'})])
    samplesDict = sample_prep(samples)
    # add 'demultiplexedPath' and 'trimmedPath' fields for each sample
    for sample in samplesDict:
        demultiplexedPath = '{experiment}/raw_data/{sampleName}.fastq.gz'.format(
            experiment=experiment,
            sampleName=samplesDict[sample]['name'])
        trimmedPath = '{experiment}/{sampleName}_trimmed.fastq'.format(
            experiment=experiment,
            sampleName=samplesDict[sample]['name'])
        samplesDict[sample]['demultiplexedPath'] = demultiplexedPath
        samplesDict[sample]['trimmedPath'] = trimmedPath

    #demultiplex based on barcodes defined in the sample file
    demultiplex_fastq(reads, samplesDict, experiment)

    # Prepare genome files from the GenBank input
    gbk2fna(gbkfile, organism, genomeDir)
    gbk2ftt(gbkfile, organism, genomeDir)

    # Change directory, build bowtie indexes, change directory back
    with cd(genomeDir):
        bowtieBuild(organism)

    # Dictionary of each sample's cpm by gene
    geneMappings = {}
    for sample in samplesDict:
        s = samplesDict[sample]
        trim_fastq(s['demultiplexedPath'], s['trimmedPath'], sample)
        # Change directory, map to bowtie, change directory back
        trimmedSampleFile = '{0}_trimmed.fastq'.format(sample)
        bowtieOutputFile = '{0}_bowtie.txt'.format(sample)
        with cd(genomeDir):
            # Paths are relative to the genome_lookup directory
            # from where bowtie is called
            bowtie_in = '../{0}'.format(trimmedSampleFile)
            bowtie_out = '../{0}'.format(bowtieOutputFile)
            # map to bowtie and produce the output file
            bowtieMap(organism, bowtie_in, bowtie_out)
        # Delete trimmed fastq file after writing mapping results
        os.remove(s['trimmedPath'])
        mapSites('{0}/{1}'.format(experiment, bowtieOutputFile))
        # Add gene-level results for the sample to geneMappings
        geneMappings[sample] = mapGenes(organism, sample, experiment)
    buildGeneTable(organism, samplesDict, geneMappings, experiment)

if __name__ == '__main__':
    main()
