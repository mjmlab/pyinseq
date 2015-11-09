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

    # List of file paths of the demultiplexed files
    demultiplexedSample_list = demultiplexedSamplesToProcess(samples, experiment)

    # Change directory, build bowtie indexes, change directory back
    with cd(genomeDir):
        bowtieBuild(organism)

    # PROCESS SAMPLE
    # Rewrite as new .fastq file with only chromosome sequence
    # Trim barcode, trim transposon. Trim corresponding quality.
    # Ignore if not a good transposon junction.
    for samplePath in demultiplexedSample_list:
        s1 = samplePath.split('.')[0].rfind('/')
        sampleName = samplePath.split('.')[0][s1+1:]
        trim_fastq(samplePath, sampleName, experiment)
        # Change directory, map to bowtie, change directory back
        with cd(genomeDir):
#            bowtieMap(organism, readsAssigned, bowtieOutput)
            pass


    ## *** RUN BOWTIE ON NEW FILE ***
    ## *** DELETE TRIMMED FASTQ FILE? YES UNLESS (K)EEP TEMPORARY FILES SELECTED ***

    ## THEN DO MATH ON EACH LINE OF THE BOWTIE OUTPUT:
    ## {(contig, position) : [Lcount, Rcount]}
    ## d.setdefault((contig,position),[0,0])[0] += 1   # L
    ## d.setdefault((contig,position),[0,0])[1] += 2   # R



    # Assign and trim barcodes
    # Now currently filtering by default (16-17 bp, TA at end)
    # In future could add as command line option
    #assignAndTrim(reads, samples, experiment, tempDir)

    # Prepare the bowtie indexes
    # Map the reads using bowtie_map
#    readsAssigned = '{0}_assigned.fastq'.format(experiment) # already exsists
#    bowtieOutput = '{0}_bowtieOutput.txt'.format(experiment) # will get created at next step


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
