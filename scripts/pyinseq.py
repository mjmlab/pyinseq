#!/usr/bin/env python
"""Main script for running the pyinseq package."""

import argparse
import os
from shutil import copyfile
import sys
import yaml
from demultiplex import sample_prep, demultiplex_fastq, trim_fastq
from gbkconvert import gbk2fna, gbk2ftt
from mapReads import bowtieBuild, bowtieMap, parseBowtie
from processMapping import mapSites, mapGenes, buildGeneTable
from utils import convert_to_filename, createExperimentDirectories

def parseArgs(args):
    """Parse command line arguments."""
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
    parser.add_argument('--nobarcodes',
                        help='barcodes have already been removed from the samples; \
                        -i should list the directory with filenames (.fastq.gz) \
                        corresponding to the sample names',
                        action='store_true',
                        default=False)
    parser.add_argument('--keepall',
                        help='keep all intermediate files generated \
                        (warning: large size!)',
                        action='store_true',
                        default=False)
    return parser.parse_args(args)


class cd:
    """Context manager to change to the specified directory then back."""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def pipeline_organize(samples):

    print('\n===================='\
          '\n*    Setting up    *'\
          '\n====================\n')

    # Create the directory struture based on the experiment name
    createExperimentDirectories(experiment)

    # Note: barcode length hardcoded at 4 bp here
    barcode_qc, barcode_length = True, 4

    # if nobarcodes:
        # barcode_qc, barcode_length = False, 0

    # TODO(For rerunning samples, modify samplesDict construction; read in a YAML file?)

    # TODO(Modify as needed for already-demultiplexed samples)

    # samples = OrderedDict([('name1', {'name': 'name1', 'barcode': 'barcode1'}),
    #    ('name2', {'name': 'name2', 'barcode': 'barcode2'})])
    global samplesDict
    samplesDict = sample_prep(samples, barcode_qc)

    # add 'demultiplexedPath' and 'trimmedPath' fields for each sample
    for sample in samplesDict:
        demultiplexedPath = 'results/{experiment}/raw_data/{sampleName}.fastq.gz'.format(
            experiment=experiment,
            sampleName=samplesDict[sample]['name'])
        trimmedPath = 'results/{experiment}/{sampleName}_trimmed.fastq'.format(
            experiment=experiment,
            sampleName=samplesDict[sample]['name'])
        samplesDict[sample]['demultiplexedPath'] = demultiplexedPath
        samplesDict[sample]['trimmedPath'] = trimmedPath

    print('\nProcessing {} total samples:'.format(len(samplesDict)))
    for s in samplesDict:
        print('{0}\n  barcode: {1}'.format(s, samplesDict[s]['barcode']))
    samples_yaml = 'results/{}/samples.yml'.format(experiment)
    with open(samples_yaml, 'w') as fo:
        fo.write(yaml.dump(samplesDict, default_flow_style=False))
    print('Sample details written to {}'.format(samples_yaml))

def pipeline_no_demultiplex(reads):
    # copy reads files into the experiment/raw_data directory
    for sample in samplesDict:
        # makes sure the reads directory has a trailing slash
        if reads[-1] != '/':
            reads += '/'
        src = reads + sample + '.fastq.gz'
        dst = samplesDict[sample]['demultiplexedPath']
        copyfile(src, dst)

def pipeline_demultiplex(reads):

    print('\n===================='\
          '\n*  Demultiplexing  *'\
          '\n====================\n')

    # demultiplex based on barcodes defined in the sample file
    print('\nDemultiplexing from input file:\n  {}'.format(reads))
    nreads = demultiplex_fastq(reads, samplesDict, experiment)
    logdata['total_reads'] = nreads
    print('Demultiplexed into output files:')
    for s in samplesDict:
        print('  ' + samplesDict[s]['demultiplexedPath'])

def pipeline_mapping(gbkfile, organism, genomeDir, disruption, barcode_length=4):
    # Prepare genome files from the GenBank input

    print('\n===================='\
          '\n*     Mapping      *'\
          '\n====================\n')

    fnaPrint = \
        '\nPreparing nucleotide fasta file from GenBank file to use in bowtie mapping.\n' \
        '  GenBank source file:  {}'.format(gbkfile)
    fttPrint = \
        '\nPreparing feature table file from GenBank file to use in gene mapping.\n' \
        '  GenBank source file:  {}'.format(gbkfile)
    print(fnaPrint)
    gbk2fna(gbkfile, organism, genomeDir)
    print(fttPrint)
    gbk2ftt(gbkfile, organism, genomeDir)

    # Change directory, build bowtie indexes, change directory back
    with cd(genomeDir):
        print('\nBuilding bowtie index files in results/{}/genome_lookup'.format(experiment))
        bowtieBuild(organism)

    # Dictionary of each sample's cpm by gene
    geneMappings = {}
    for sample in samplesDict:
        s = samplesDict[sample]
        print('\nProcessing sample {}'.format(sample))
        sample_reads, trimmed_reads = trim_fastq(s['demultiplexedPath'], s['trimmedPath'], sample, barcode_length)
        logdata[sample] = {}
        logdata[sample]['reads_with_bc'] = sample_reads
        logdata[sample]['reads_with_bc_seq_tn'] = trimmed_reads
        # Change directory, map to bowtie, change directory back
        trimmedSampleFile = '{0}_trimmed.fastq'.format(sample)
        bowtieOutputFile = '{0}_bowtie.txt'.format(sample)
        with cd(genomeDir):
            # Paths are relative to the genome_lookup directory
            # from where bowtie is called
            bowtie_in = '../{0}'.format(trimmedSampleFile)
            bowtie_out = '../{0}'.format(bowtieOutputFile)
            # map to bowtie and produce the output file
            print('\nMapping {} reads with bowtie'.format(sample))
            bowtie_msg_out = bowtieMap(organism, bowtie_in, bowtie_out)
            # store bowtie data for each sample in dictionary
            logdata[sample]['bowtie_results'] = parseBowtie(bowtie_msg_out)
        # Map each bowtie result to the chromosome
        insertions = len(mapSites('results/{0}/{1}'.format(experiment, bowtieOutputFile)))
        logdata[sample]['insertion_sites'] = insertions
        # Add gene-level results for the sample to geneMappings
        # Filtered on gene fraction disrupted as specified by -d flag
        geneMappings[sample] = mapGenes(organism, sample, disruption, experiment)
        if not keepall:
            # Delete trimmed fastq file, bowtie mapping file after writing mapping results
            os.remove(s['trimmedPath'])
            os.remove('results/{0}/{1}'.format(experiment, bowtieOutputFile))
    buildGeneTable(organism, samplesDict, geneMappings, experiment)
    # print(logdata)


def pipeline_analysis():

    print('\n===================='\
          '\n*     Analysis     *'\
          '\n====================\n')

    samples_summary = 'results/{}/samples_summary.yml'.format(experiment)
    with open(samples_summary, 'w') as fo:
        fo.write(yaml.dump(logdata, default_flow_style=False))
    print('Writing file with summary of results:\n  {}'.format(samples_summary))

def main():
    """Start here."""
    args = parseArgs(sys.argv[1:])
    global experiment
    experiment = convert_to_filename(args.experiment)
    gbkfile = args.genome
    reads = args.input
    samples = args.samples
    # TODO(Test that disruption is between 0.0 and 1.0 (or absent, default 1.0))
    disruption = float(args.disruption)
    nobarcodes = args.nobarcodes
    global keepall
    keepall = args.keepall
    # Logging of sample info
    global logdata
    logdata = {}
    # Organism reference files called 'genome.fna' etc
    organism = 'genome'

    # --- ORGANIZE SAMPLE LIST AND FILE PATHS --- #
    pipeline_organize(samples)

    # --- DEMULTIPLEX OR MOVE FILES IF ALREADY DEMULTIPLEXED --- #
    if nobarcodes:
        pipeline_no_demultiplex(reads)
    else:
        pipeline_demultiplex(reads)

    # --- BOWTIE MAPPING --- #
    genomeDir = 'results/{experiment}/genome_lookup/'.format(experiment=experiment)
    pipeline_mapping(gbkfile, organism, genomeDir, disruption)

    # --- ANALYSIS OF RESULTS --- #
    pipeline_analysis()


    # --- CONFIRM COMPLETION --- #
    print('\n===================='\
          '\n*       Done       *'\
          '\n====================\n')


if __name__ == '__main__':
    main()
