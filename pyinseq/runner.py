#!/usr/bin/env python3

"""Main script for running the pyinseq package."""
import argparse
import csv
import glob
import logging
import os
import numpy as np
import pandas as pd
import re
import screed
import subprocess
import sys
import yaml
from shutil import copyfile
from collections import OrderedDict
from .analyze import read_sites_file, nfifty, plot_insertions
from .demultiplex import demultiplex_fastq, write_reads
from .gbkconvert import gbk2fna, gbk2ftt
from .mapReads import bowtie_build, bowtie_map, parse_bowtie
from .processMapping import map_sites, map_genes, build_gene_table
from .utils import (
    convert_to_filename,
    create_experiment_directories,
)  # has logging config

# Note: stdout logging is set in utils.py
logger = logging.getLogger("pyinseq")
logger.setLevel(logging.INFO)


def parseArgs(args):
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", help="input Illumina reads file or folder", required=True
    )
    parser.add_argument(
        "-s", "--samples", help="sample list with barcodes", required=True
    )
    parser.add_argument(
        "-e",
        "--experiment",
        help="experiment name (no spaces or special characters)",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--genome",
        help="genome in GenBank format (one concatenated file for multiple contigs/chromosomes)",
        required=True,
    )
    parser.add_argument(
        "-d", "--disruption", help="fraction of gene disrupted (0.0 - 1.0)", default=1.0
    )
    parser.add_argument(
        "--min_count", help="Minimum number of reads per insertion site", default=3
    )
    parser.add_argument(
        "--max_ratio",
        help="Maximum ratio of left:right or right:left reads per insertion site",
        default=10,
    )
    """Inactive arguments in current version
    parser.add_argument('-s', '--samples',
                        help='sample list with barcodes. \
                        If not provided then entire folder provided for --input is analyzed',
                        required=False)
    parser.add_argument('--nobarcodes',
                        help='barcodes have already been removed from the samples; \
                        -i should list the directory with filenames (.fastq.gz) \
                        corresponding to the sample names',
                        action='store_true',
                        default=False)
    parser.add_argument('--compress',
                        help='compress (gzip) demultiplexed samples',
                        action='store_true',
                        default=False)
    parser.add_argument('--keepall',
                        help='keep all intermediate files generated',
                        action='store_true',
                        default=False)"""
    return parser.parse_args(args)


def demultiplex_parseArgs(args):
    """Parse command line arguments for `pyinseq demultiplex`."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", help="input Illumina reads file or folder", required=True
    )
    parser.add_argument(
        "-s", "--samples", help="sample list with barcodes", required=True
    )
    parser.add_argument(
        "-e",
        "--experiment",
        help="experiment name (no spaces or special characters)",
        required=True,
    )
    parser.add_argument(
        "--notrim",
        help="do not write trimmed reads (i.e. write raw reads only)",
        action="store_true",
        required=False,
    )
    return parser.parse_args(args)


def genomeprep_parseArgs(args):
    """Parse command line arguments for `pyinseq genomeprep`."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-e",
        "--experiment",
        help="experiment name (no spaces or special characters)",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--genome",
        help="genome in GenBank format (one concatenated file for multiple contigs/chromosomes)",
        required=True,
    )
    parser.add_argument(
        "--noindex",
        help="do not generate bowtie indexes",
        action="store_true",
        required=False,
    )
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


class Settings:
    """Instantiate to set up settings for the experiment."""

    def __init__(self, experiment_name):
        # command options are: ['pyinseq', 'demultiplex']
        self.command = "pyinseq"
        self.experiment = convert_to_filename(experiment_name)
        self.path = "results/{}/".format(self.experiment)
        self.parse_genbank_file = True
        self.genome_path = self.path + "genome_lookup/"
        # organism reference files called 'genome.fna' etc
        self.organism = "genome"
        self.raw_path = self.path + "raw_data/"
        self.generate_bowtie_index = True
        self.process_reads = True
        self.process_sample_list = True
        self.map_to_genome = True
        self.samples_yaml = self.path + "samples.yml"
        self.summary_log = self.path + "log.txt"
        self.write_trimmed_reads = True
        self.min_counts = 3  # counts at one transposon site for it to qualify
        self.max_ratio = 10  # max ratio of left/right sites for it to qualify
        # may be modified
        self.keepall = False
        self.barcode_length = 4
        self.disruption = 1

    def __repr__(self):
        # Print each variable on a separate line
        return "\n".join("%s: %s" % item for item in vars(self).items())

    def set_command_specific_settings(self, cmd):
        if cmd == "demultiplex":
            self.command = "demultiplex"
            self.parse_genbank_file = False
            self.generate_bowtie_index = False
            self.map_to_genome = False
        if cmd == "genomeprep":
            self.command = "genomeprep"
            self.process_reads = False
            self.process_sample_list = False
            self.map_to_genome = False


def set_paths(experiment_name):
    """Set up relative paths for subdirectories."""
    experiment = convert_to_filename(experiment_name)
    samples_yaml = "results/{}/samples.yml".format(experiment)
    summary_log = "results/{}/log.txt".format(experiment)
    path = {
        "experiment": experiment,
        "samples_yaml": samples_yaml,
        "summary_log": summary_log,
    }
    return path


def set_disruption(d, setting):
    """Check that gene disruption is 0.0 to 1.0; otherwise set to 1.0."""
    if d < 0.0 or d > 1.0:
        logger.error(
            "Disruption value provided ({0}) is not in range 0.0 to 1.0; proceeding with default value of 1.0".format(
                d
            )
        )
        d = 1.0

    setting.disruption = d
    return


def set_gene_parameters(min_count, max_ratio, setting):
    """Check that min_count and max_ratio are positive, otherwise set default."""
    if min_count < 0:
        logger.error(
            "Min_count value provided ({0}) is not positive; proceeding with default value of 3".format(
                min_count
            )
        )
    else:
        setting.min_counts = min_count

    if max_ratio < 0:
        logger.error(
            "Max_ratio value provided ({0}) is not positive; proceeding with default value of 10".format(
                max_ratio
            )
        )
    else:
        setting.max_ratio = max_ratio
    return


def tab_delimited_samples_to_dict(sample_file):
    """Read sample names, barcodes from tab-delimited into an OrderedDict."""
    samplesDict = OrderedDict()
    with open(sample_file, "r", newline="") as csvfile:
        for line in csv.reader(csvfile, delimiter="\t"):
            # ignore comment lines in original file
            if not line[0].startswith("#"):
                # sample > filename-acceptable string
                # barcode > uppercase
                sample = convert_to_filename(line[0])
                barcode = line[1].upper()
                if sample not in samplesDict and barcode not in samplesDict.values():
                    samplesDict[sample] = {"barcode": barcode}
                else:
                    raise IOError(
                        "Error: duplicate sample {0} barcode {1}".format(
                            sample, barcode
                        )
                    )
    return samplesDict


def yaml_samples_to_dict(sample_file):
    """Read sample names, barcodes from yaml into an OrderedDict."""
    with open(sample_file, "r") as f:
        samplesDict = yaml.load(sample_file)
    return samplesDict


def directory_of_samples_to_dict(directory):
    """Read sample names from a directory of .gz files into an OrderedDict."""
    samplesDict = OrderedDict()
    for gzfile in list_files(directory):
        # TODO(convert internal periods to underscore? use regex?)
        # extract file name before any periods
        f = os.path.splitext(os.path.basename(gzfile))[0].split(".")[0]
        samplesDict[f] = {}
    return samplesDict


def list_files(folder, ext="gz"):
    """Return list of .gz files from the specified folder."""
    with cd(folder):
        return [f for f in glob.glob("*.{}".format(ext))]


def build_fna_and_ftt_files(gbkfile, settings):
    """Convert GenBank file to a fasta nucleotide and feature table files."""
    gbk2fna(gbkfile, settings.organism, settings.genome_path)
    gbk2ftt(gbkfile, settings.organism, settings.genome_path)


def build_bowtie_index(settings):
    """Change directory, build bowtie indexes, and change directory back."""
    with cd(settings.genome_path):
        logger.info(
            "Building bowtie index files in results/{}/genome_lookup".format(
                settings.experiment
            )
        )
        bowtie_build(settings.organism)


def pipeline_mapping(settings, samplesDict):
    """Aggregate bowtie output, map to genes in the feature table, and aggregate samples."""
    # Dictionary of each sample's cpm by gene
    geneMappings = {}
    mapping_data = {}
    for sample in samplesDict:
        with cd(settings.genome_path):
            # Paths are relative to the genome_lookup directory
            # from where bowtie is called
            bowtie_in = "../" + sample + "_trimmed.fastq"
            bowtie_out = "../" + sample + "_bowtie.txt"
            # map to bowtie and produce the output file
            logger.info("Sample {}: map reads with bowtie".format(sample))
            bowtie_msg_out = bowtie_map(settings.organism, bowtie_in, bowtie_out)
            logger.info(bowtie_msg_out)
            # store bowtie data for each sample in dictionary
            mapping_data[sample] = {"bowtie_results": [], "insertion_sites": []}
            mapping_data[sample]["bowtie_results"] = parse_bowtie(bowtie_msg_out)
        # Map each bowtie result to the chromosome
        logger.info(
            "Sample {}: summarize the site data from the bowtie results".format(sample)
        )
        insertions = len(map_sites(sample, samplesDict, settings))
        mapping_data[sample]["insertion_sites"] = insertions
        # Add gene-level results for the sample to geneMappings
        # Filtered on gene fraction disrupted as specified by -d flag
        logger.info("Sample {}: map site data to genes".format(sample))
        geneMappings[sample] = map_genes(sample, settings)
        # if not settings.keepall:
        #    # Delete trimmed fastq file, bowtie mapping file after writing mapping results
        #    os.remove(s['trimmedPath'])
        #    os.remove('results/{0}/{1}'.format(Settings.experiment, bowtieOutputFile))
        nfifty_result = nfifty(sample, settings)
        logger.info("N50 result for {}: {}".format(sample, nfifty_result))
    logger.info("Aggregate gene mapping from all samples into the summary_data_table")
    build_gene_table(settings.organism, samplesDict, geneMappings, settings.experiment)


def pipeline_summarize(samplesDict, settings, typed_command_after_pyinseq):
    """Summary of INSeq run."""
    logger.info("Print samples info: {}".format(settings.samples_yaml))
    with open(settings.samples_yaml, "w") as fo:
        fo.write(yaml.dump(samplesDict, default_flow_style=False))

    # write summary log with more data
    logger.info("Print summary log: {}".format(settings.summary_log))
    logger.info(
        "Print command entered" + "\npyinseq " + " ".join(typed_command_after_pyinseq)
    )
    logger.info("Print settings" + "\n" + str(settings))
    logger.info(
        "Print samples detail" + "\n" + yaml.dump(samplesDict, default_flow_style=False)
    )


# def pipeline_analysis(samplesDict, settings):
#    """Analysis of output."""
#    for sample in samplesDict:
#        print('N50', sample, nfifty(sample, settings))
#        # plot_insertions(sample, settings)


def main(args):
    """Start here."""
    logger.info("Process command line arguments")
    typed_command_after_pyinseq = args
    # pyinseq with nothing typed after it
    if args == []:
        args = ["-h"]
        # command = 'pyinseq'
        # args = parseArgs(args)
    # pyinseq demultiplex
    if args[0] == "demultiplex":
        command = "demultiplex"
        args = demultiplex_parseArgs(args[1:])
    elif args[0] == "genomeprep":
        command = "genomeprep"
        args = genomeprep_parseArgs(args[1:])
    else:
        command = "pyinseq"
        args = parseArgs(args)
    # Initialize the settings object
    settings = Settings(args.experiment)
    settings.set_command_specific_settings(command)
    try:
        # for `pyinseq demultiplex` only
        settings.write_trimmed_reads = not args.notrim
    except:
        pass
    try:
        # for `pyinseq genomeprep` only
        settings.generate_bowtie_index = not args.noindex
        print("args.noindex", args.noindex)
        print("settings.generate_bowtie_index", settings.generate_bowtie_index)
    except:
        pass
    # Keep intermediate files
    settings.keepall = False  # args.keepall
    if settings.process_reads:
        reads = args.input
    if settings.parse_genbank_file:
        gbkfile = args.genome
        if settings.process_reads:
            set_disruption(float(args.disruption), settings)
            set_gene_parameters(int(args.min_count), int(args.max_ratio), settings)
    # sample names and paths
    if settings.process_sample_list:
        samples = args.samples
        # barcodes_present = not args.nobarcodes
        if samples:
            samplesDict = tab_delimited_samples_to_dict(samples)
        else:
            reads = os.path.abspath(reads)
            samplesDict = directory_of_samples_to_dict(samples)
    try:
        logger.debug("samplesDict: {0}".format(samplesDict))
    except (UnboundLocalError):
        pass

    # --- SET UP DIRECTORIES --- #
    create_experiment_directories(settings)

    # --- SET UP LOG FILE --- #
    fh = logging.FileHandler(settings.summary_log)
    fh.setFormatter(
        logging.Formatter(
            "%(asctime)s - %(levelname)s - %(module)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M",
        )
    )  # also set in utils
    logger.addHandler(fh)

    # --- WRITE DEMULTIPLEXED AND TRIMMED FASTQ FILES --- #
    if settings.process_reads:
        logger.info("Demultiplex reads")
        demultiplex_fastq(reads, samplesDict, settings)

    # --- MAPPING TO SITES AND GENES --- #
    if settings.parse_genbank_file:
        logger.info("Prepare genome features (.ftt) and fasta nucleotide (.fna) files")
        build_fna_and_ftt_files(gbkfile, settings)
    if settings.generate_bowtie_index:
        logger.info("Prepare bowtie index")
        build_bowtie_index(settings)
    if settings.map_to_genome:
        logger.info("Map with bowtie")
        pipeline_mapping(settings, samplesDict)

    # if not samples:
    #    Settings.summaryDict['total reads'] = 0
    #    for sample in Settings.samplesDict:
    #        print(Settings.samplesDict[sample])
    #        Settings.summaryDict['total reads'] += Settings.samplesDict[sample]['reads_with_bc']

    # --- SUMMARY OF RESULTS --- #
    if settings.command in ["pyinseq"]:
        pipeline_summarize(samplesDict, settings, typed_command_after_pyinseq)

    # --- CONFIRM COMPLETION --- #
    logger.info("***** {} complete! *****".format(settings.command))


if __name__ == "__main__":
    main()
