#!/usr/bin/env python3

'''Main script for running the pyinseq package.'''
import argparse
import csv
import glob
import logging
import os
import numpy as np
import pandas as pd
import regex as re
import screed
import sys
import yaml
from shutil import copyfile
from collections import OrderedDict
from .analyze import read_sites_file, nfifty, plot_insertions
from .demultiplex import demultiplex_fastq, write_reads
from .gbkconvert import gbk2fna, gbk2ftt
from .mapReads import bowtie_build, bowtie_map, parse_bowtie
from .processMapping import map_sites, map_genes, build_gene_table
from .utils import convert_to_filename, create_experiment_directories # has logging config

logger = logging.getLogger(__name__)

def parseArgs(args):
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',
                        help='input Illumina reads file or folder',
                        required=True)
    parser.add_argument('-s', '--samples',
                        help='sample list with barcodes. \
                        If not provided then entire folder provided for --input is analyzed',
                        required=False)
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
    parser.add_argument('--demultiplex',
                        help='demultiplex initial file into separate files by barcode',
                        action='store_true',
                        default=False)
    parser.add_argument('--compress',
                        help='compress (gzip) demultiplexed samples',
                        action='store_true',
                        default=False)
    parser.add_argument('--keepall',
                        help='keep all intermediate files generated',
                        action='store_true',
                        default=False)
    return parser.parse_args(args)


class cd:
    '''Context manager to change to the specified directory then back.'''
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

class Settings():
    '''Instantiate to set up settings for the experiment'''
    def __init__(self, experiment_name):
        # standard
        self.experiment = convert_to_filename(experiment_name)
        self.path = 'results/{}/'.format(self.experiment)
        self.genome_path = self.path + 'genome_lookup/'
        self.genome_index_path = self.path + 'genome_lookup/genome'
        self.raw_path = self.path + 'raw_data/'
        self.samples_yaml = self.path + 'samples.yml'
        self.summary_yaml = self.path + 'summary.yml'
        # may be modified
        self.keepall = False
        self.barcode_length = 4

    # Set up directories?


def set_paths(experiment_name):
    experiment = convert_to_filename(experiment_name)
    samples_yaml = 'results/{}/samples.yml'.format(experiment)
    summary_yaml = 'results/{}/summary.yml'.format(experiment)
    path = {'experiment': experiment,
            'samples_yaml': samples_yaml,
            'summary_yaml': summary_yaml}
    return path


def set_disruption(d):
    '''Check that gene disrution is 0.0 to 1.0; otherwise set to 1.0'''
    if d < 0.0 or d > 1.0:
        logger.error('Disruption value provided ({0}) is not in range 0.0 to 1.0; proceeding with default value of 1.0'.format(d))
        d = 1.0
    return d


def tab_delimited_samples_to_dict(sample_file):
    '''Read sample names, barcodes from tab-delimited into an OrderedDict.'''
    samplesDict = OrderedDict()
    with open(sample_file, 'r', newline='') as csvfile:
        for line in csv.reader(csvfile, delimiter='\t'):
            if not line[0].startswith('#'):  # ignore comment lines in original file
                # sample > filename-acceptable string
                # barcode > uppercase
                sample = convert_to_filename(line[0])
                barcode = line[1].upper()
                if sample not in samplesDict and barcode not in samplesDict.values():
                        samplesDict[sample] = {'barcode': barcode}
                else:
                    raise IOError('Error: duplicate sample {0} barcode {1}'.format(sample, barcode))
    return samplesDict


def yaml_samples_to_dict(sample_file):
    '''Read sample names, barcodes from yaml into an OrderedDict.'''
    with open(sample_file, 'r') as f:
        samplesDict = yaml.load(sample_file)
    return samplesDict


def directory_of_samples_to_dict(directory):
    '''Read sample names from a directory of .gz files into an OrderedDict.'''
    samplesDict = OrderedDict()
    for gzfile in list_files(directory):
        # TODO(convert internal periods to underscore? use regex?)
        # extract file name before any periods
        f = (os.path.splitext(os.path.basename(gzfile))[0].split('.')[0])
        samplesDict[f] = {}
    return samplesDict


def list_files(folder, ext='gz'):
    '''Return list of .gz files from the specified folder'''
    with cd(folder):
        return [f for f in glob.glob('*.{}'.format(ext))]


def build_fna_and_ftt_files(gbkfile, organism, settings):
    gbk2fna(gbkfile, organism, settings.genome_path)
    gbk2ftt(gbkfile, organism, settings.genome_path)


def build_bowtie_index(organism, settings):
    # Change directory, build bowtie indexes, change directory back
    with cd(settings.genome_path):
        logger.info('Building bowtie index files in results/{}/genome_lookup'.format(settings.experiment))
        bowtie_build(organism)


def pipeline_mapping(organism, settings, samplesDict, disruption):
    # Dictionary of each sample's cpm by gene
    geneMappings = {}
    mapping_data = {}
    for sample in samplesDict:
        with cd(settings.genome_path):
            # Paths are relative to the genome_lookup directory
            # from where bowtie is called
            bowtie_in = '../' + sample + '_trimmed.fastq'
            bowtie_out = '../' + sample + '_bowtie.txt'
            # map to bowtie and produce the output file
            logger.info('Sample {}: map reads with bowtie'.format(sample))
            bowtie_msg_out = bowtie_map(organism, bowtie_in, bowtie_out)
            # store bowtie data for each sample in dictionary
            mapping_data[sample] = {'bowtie_results': [], 'insertion_sites': []}
            mapping_data[sample]['bowtie_results'] = parse_bowtie(bowtie_msg_out)
        # Map each bowtie result to the chromosome
        logger.info('Sample {}: summarize the site data from the bowtie results'.format(sample))
        insertions = len(map_sites(sample, samplesDict, settings))
        mapping_data[sample]['insertion_sites'] = insertions
        # Add gene-level results for the sample to geneMappings
        # Filtered on gene fraction disrupted as specified by -d flag
        logger.info('Sample {}: map site data to genes'.format(sample))
        geneMappings[sample] = map_genes(organism, sample, disruption, settings)
        #if not settings.keepall:
        #    # Delete trimmed fastq file, bowtie mapping file after writing mapping results
        #    os.remove(s['trimmedPath'])
        #    os.remove('results/{0}/{1}'.format(Settings.experiment, bowtieOutputFile))
    logger.info('Aggregate gene mapping from all samples into the summary_data_table'.format(sample))
    build_gene_table(organism, samplesDict, geneMappings, settings.experiment)


def pipeline_analysis(samplesDict, settings):
    # logger.info('Print summary logs.')
    # print('Writing file with summary data for each sample:\n  {}'.format(settings.samples_yaml))
    # print(settings.samples_yaml)
    # logger.info('samplesDict: {}'.format(samplesDict))
    # print(yaml.dump(settings.samples_yaml, default_flow_style=False))
    # with open(settings.samples_yaml, 'w') as fo:
    #    fo.write(yaml.dump(samplesDict, default_flow_style=False))

    # write summary.yml with more data
    # print('Writing file with overall summary information:\n  {}'.format(settings.summary_yaml))
    # print(yaml.dump(settings.summary_yaml, default_flow_style=False))
    # with open(settings.summary_yaml, 'w') as fo:
    #    fo.write(yaml.dump(settings.summary_yaml, default_flow_style=False))
    # analyze individual samples
    for sample in samplesDict:
        print('N50', sample, nfifty(sample, settings))
        # plot_insertions(sample, settings)


def main(args):
    '''Start here.'''
    logger.info('Process command line arguments')
    args = parseArgs(args)
    # Initialize the settings object
    settings = Settings(args.experiment)
    # Keep intermediate files
    settings.keepall = args.keepall
    gbkfile = args.genome
    reads = args.input
    disruption = set_disruption(float(args.disruption))
    # Organism reference files called 'genome.fna' etc
    organism = 'genome'
    # sample names and paths
    samples = args.samples
    barcodes_present = not args.nobarcodes
    if samples:
        samplesDict = tab_delimited_samples_to_dict(samples)
    else:
        reads = os.path.abspath(reads)
        samplesDict = directory_of_samples_to_dict(samples)
    logger.debug('samplesDict: {0}'.format(samplesDict))

    # --- SET UP DIRECTORIES --- #
    create_experiment_directories(settings.experiment)

    # --- WRITE DEMULTIPLEXED AND TRIMED FASTQ FILES --- #
    logger.info('Demultiplex reads')
    demultiplex_fastq(reads, samplesDict, settings)

    # --- MAPPING TO SITES AND GENES --- #
    logger.info('Prepare genome features (.ftt) and fasta nucleotide (.fna) files')
    build_fna_and_ftt_files(gbkfile, organism, settings)
    logger.info('Prepare bowtie index')
    build_bowtie_index(organism, settings)
    logger.info('Map with bowtie')
    pipeline_mapping(organism, settings, samplesDict, disruption)

    #if not samples:
    #    Settings.summaryDict['total reads'] = 0
    #    for sample in Settings.samplesDict:
    #        print(Settings.samplesDict[sample])
    #        Settings.summaryDict['total reads'] += Settings.samplesDict[sample]['reads_with_bc']

    # --- ANALYSIS OF RESULTS --- #
    pipeline_analysis(samplesDict, settings)

    # --- CONFIRM COMPLETION --- #
    logger.info('***** pyinseq pipeline complete! *****')

if __name__ == '__main__':
    main()
