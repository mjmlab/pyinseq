#!/usr/bin/env python3

"""Main script for running the pyinseq package."""
import argparse
import csv
import glob
import logging
import os
import yaml

from .analyze import t_fifty, spearman_correlation
from .demultiplex import demultiplex_fastq
from .gbk_convert import gbk2fna, gbk2ftt
from .map_reads import bowtie_build, bowtie_map, parse_bowtie
from .process_mapping import map_sites, map_genes, build_gene_table
from .utils import (
    convert_to_filename,
    create_experiment_directories,
)  # has logging config

# Note: stdout logging is set in utils.py
logger = logging.getLogger("pyinseq")


def parse_args(args):
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
    parser.add_argument('--keep_all',
                        help='keep all intermediate files generated',
                        action='store_true',
                        default=False)"""
    return parser.parse_args(args)


def demultiplex_parse_args(args):
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


def genomeprep_parse_args(args):
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

    def __init__(self, new_path):
        self.new_path = os.path.expanduser(new_path)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.new_path)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


class Settings:
    """Instantiate to set up settings for the experiment."""

    def __init__(self, experiment_name):
        # command options are: ['pyinseq', 'demultiplex']
        self.command = "pyinseq"
        self.experiment = convert_to_filename(experiment_name)
        self.path = f"results/{self.experiment}/"
        self.parse_genbank_file = True
        self.genome_path = f"{self.path}genome_lookup/"
        # organism reference files called 'genome.fna' etc
        self.organism = "genome"
        self.raw_path = f"{self.path}raw_data/"
        # self.figures_path = f"{self.path}figures/"
        # self.analysis_path = f"{self.path}analysis/"
        self.generate_bowtie_index = True
        self.process_reads = True
        self.process_sample_list = True
        self.map_to_genome = True
        self.samples_yaml = f"{self.path}samples.yml"
        self.summary_log = f"{self.path}log.txt"
        self.write_trimmed_reads = True
        self.summary_table = f"{self.path}summary_gene_table.txt"
        self.min_counts = 3  # counts at one transposon site for it to qualify
        self.max_ratio = 10  # max ratio of left/right sites for it to qualify
        # may be modified
        self.keep_all = False
        self.barcode_length = 4
        self.disruption = 1

    def __repr__(self):
        # Print each variable on a separate line
        return "\n".join("%s: %s" % item for item in vars(self).items())

    def _set_command_specific_settings(self, cmd):
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


def _set_paths(experiment_name):
    """Set up relative paths for subdirectories."""
    experiment = convert_to_filename(experiment_name)
    samples_yaml = f"results/{experiment}/samples.yml"
    summary_log = f"results/{experiment}/log.txt"
    path = {
        "experiment": experiment,
        "samples_yaml": samples_yaml,
        "summary_log": summary_log,
    }
    return path


def _set_disruption(d, setting):
    """Check that gene disruption is 0.0 to 1.0; otherwise set to 1.0."""
    if d < 0.0 or d > 1.0:
        logger.error(
            f"Disruption value provided ({d}) is not in range 0.0 to 1.0; proceeding with default value of 1.0"
        )
        d = 1.0

    setting.disruption = d
    return


def _set_gene_parameters(min_count, max_ratio, setting):
    """Check that min_count and max_ratio are positive, otherwise set default."""
    if min_count < 0:
        logger.error(
            f"Min_count value provided ({min_count}) is not positive; proceeding with default value of 3"
        )
    else:
        setting.min_counts = min_count

    if max_ratio < 0:
        logger.error(
            f"Max_ratio value provided ({max_ratio}) is not positive; proceeding with default value of 10"
        )
    else:
        setting.max_ratio = max_ratio
    return


def tab_delimited_samples_to_dict(sample_file):
    """Read sample names, barcodes from tab-delimited into an OrderedDict."""
    samples_dict = {}
    with open(sample_file, "r", newline="") as csv_file:
        for line in csv.reader(csv_file, delimiter="\t"):
            # ignore comment lines in original file
            if not line[0].startswith("#"):
                # sample > filename-acceptable string
                # barcode > uppercase
                sample = convert_to_filename(line[0])
                barcode = line[1].upper()
                if sample not in samples_dict and barcode not in samples_dict.values():
                    samples_dict[sample] = {"barcode": barcode}
                else:
                    raise IOError(f"Error: duplicate sample {sample} barcode {barcode}")
    return samples_dict


def yaml_samples_to_dict(sample_file):
    """Read sample names, barcodes from yaml."""
    with open(sample_file, "r") as f:
        samples_dict = yaml.load(sample_file)
    return samples_dict


def directory_of_samples_to_dict(directory):
    """Read sample names from a directory of .gz files."""
    samples_dict = {}
    for gz_file in list_files(directory):
        # TODO(convert internal periods to underscore? use regex?)
        # extract file name before any periods
        f = os.path.splitext(os.path.basename(gz_file))[0].split(".")[0]
        samples_dict[f] = {}
    return samples_dict


def list_files(folder, ext="gz"):
    """Return list of .gz files from the specified folder."""
    with cd(folder):
        return [f for f in glob.glob(f"*.{ext}")]


def build_fna_and_ftt_files(gbk_file, settings):
    """Convert GenBank file to a fasta nucleotide and feature table files."""
    gbk2fna(gbk_file, settings.organism, settings.genome_path)
    gbk2ftt(gbk_file, settings.organism, settings.genome_path)


def build_bowtie_index(settings):
    """Change directory, build bowtie indexes, and change directory back."""
    with cd(settings.genome_path):
        logger.info(
            f"Building bowtie index files in results/{settings.experiment}/genome_lookup"
        )
        bowtie_build(settings.organism)


def pipeline_mapping(settings, samples_dict):
    """Aggregate bowtie output, map to genes in the feature table, and aggregate samples."""
    # Dictionary of each sample's cpm by gene
    gene_mappings = {}
    mapping_data = {}
    for sample in samples_dict:
        with cd(settings.genome_path):
            # Paths are relative to the genome_lookup directory
            # from where bowtie is called
            bowtie_in = "../" + sample + "_trimmed.fastq"
            bowtie_out = "../" + sample + "_bowtie.txt"
            # map to bowtie and produce the output file
            logger.info(f"Sample {sample}: map reads with bowtie")
            bowtie_msg_out = bowtie_map(settings.organism, bowtie_in, bowtie_out)
            logger.info(bowtie_msg_out)
            # store bowtie data for each sample in dictionary
            mapping_data[sample] = {"bowtie_results": [], "insertion_sites": []}
            mapping_data[sample]["bowtie_results"] = parse_bowtie(bowtie_msg_out)
        # Map each bowtie result to the chromosome
        logger.info(f"Sample {sample}: summarize the site data from the bowtie results")
        insertions = len(map_sites(sample, samples_dict, settings))
        mapping_data[sample]["insertion_sites"] = insertions
        # Add gene-level results for the sample to geneMappings
        # Filtered on gene fraction disrupted as specified by -d flag
        logger.info(f"Sample {sample}: map site data to genes")
        gene_mappings[sample] = map_genes(sample, settings)
        # if not settings.keep_all:
        #    # Delete trimmed fastq file, bowtie mapping file after writing mapping results
        #    os.remove(s["trimmedPath"])
        #    os.remove("results/{Settings.experiment}/{bowtieOutputFile}"
        t_fifty_result = t_fifty(sample, settings)
        logger.info(f"T50 result for {sample}: {t_fifty_result}")
    logger.info("Aggregate gene mapping from all samples into the summary_data_table")
    build_gene_table(
        settings.organism, samples_dict, gene_mappings, settings.experiment
    )


def pipeline_summarize(samples_dict, settings, typed_command_after_pyinseq):
    """Summary of INSeq run."""
    # Yaml dump
    logger.info(f"Print samples info: {settings.samples_yaml}")
    with open(settings.samples_yaml, "w") as fo:
        fo.write(yaml.dump(samples_dict, default_flow_style=False))

    # Write summary log with more data
    logger.info(f"Print summary log: {settings.summary_log}")
    logger.info(
        "Print command entered" + "\npyinseq " + " ".join(typed_command_after_pyinseq)
    )
    logger.info("Print settings" + "\n" + str(settings))
    logger.info(
        "Print samples detail"
        + "\n"
        + yaml.dump(samples_dict, default_flow_style=False)
    )


def pipeline_analysis(samples_dict: dict, settings: Settings) -> None:
    """Analysis of output."""
    # TODO: add looping over samples and function calls from analyze script
    # T50 calculation
    T50_dict = dict()
    for sample in samples_dict:
        res_T50 = t_fifty(sample, settings)
        logger.info(f"T50 {sample}: {res_T50}")
        T50_dict[sample] = res_T50
        # plot_insertions(sample, settings)
    return


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
        args = demultiplex_parse_args(args[1:])
    # pyinseq genomeprep
    elif args[0] == "genomeprep":
        command = "genomeprep"
        args = genomeprep_parse_args(args[1:])
    # pyinseq
    else:
        command = "pyinseq"
        args = parse_args(args)
    # Initialize the settings object
    settings = Settings(args.experiment)
    settings._set_command_specific_settings(command)
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
    settings.keep_all = False  # args.keep_all
    if settings.process_reads:
        reads = args.input
    if settings.parse_genbank_file:
        gbk_file = args.genome
        if settings.process_reads:
            _set_disruption(float(args.disruption), settings)
            _set_gene_parameters(int(args.min_count), int(args.max_ratio), settings)
    # sample names and paths
    if settings.process_sample_list:
        samples = args.samples
        # barcodes_present = not args.nobarcodes
        if samples:
            samples_dict = tab_delimited_samples_to_dict(samples)
        else:
            reads = os.path.abspath(reads)
            samples_dict = directory_of_samples_to_dict(samples)
    try:
        logger.debug(f"samples_dict: {samples_dict}")
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
        demultiplex_fastq(reads, samples_dict, settings)

    # --- MAPPING TO SITES AND GENES --- #
    if settings.parse_genbank_file:
        logger.info("Prepare genome features (.ftt) and fasta nucleotide (.fna) files")
        build_fna_and_ftt_files(gbk_file, settings)
    if settings.generate_bowtie_index:
        logger.info("Prepare bowtie index")
        build_bowtie_index(settings)
    if settings.map_to_genome:
        logger.info("Map with bowtie")
        pipeline_mapping(settings, samples_dict)

    # if not samples:
    #    Settings.summaryDict['total reads'] = 0
    #    for sample in Settings.samples_dict:
    #        print(Settings.samples_dict[sample])
    #        Settings.summaryDict['total reads'] += Settings.samples_dict[sample]['reads_with_bc']

    # --- ANALYSIS OF RESULTS --- #
    if settings.command in ["pyinseq"]:
        pipeline_analysis(samples_dict, settings)

    # --- SUMMARY OF RESULTS --- #
    if settings.command in ["pyinseq"]:
        pipeline_summarize(samples_dict, settings, typed_command_after_pyinseq)

    # --- CONFIRM COMPLETION --- #
    logger.info(f"***** {settings.command} complete! *****")


if __name__ == "__main__":
    main()
