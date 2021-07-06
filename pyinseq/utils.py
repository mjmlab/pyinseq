#!/usr/bin/env python3

"""

 Functions that have no home...

"""

import re
import os
import csv
import yaml
import glob
import json
import shutil
import pickle
import subprocess as sub
from pathlib import Path

# Module imports
from pyinseq.logger import pyinseq_logger

logger = pyinseq_logger.logger


def create_experiment_directories(settings):
    """
    Create the project directory and subdirectories

    Attempt to create the directory structure:

    results/
    |
    +-{settings.experiment}/        # User-specified experiment name
      |
      +-raw_data/          # For demultiplexed reads
      |
      +-genome_lookup/     # Genome fna and ftt files, bowtie indexes

    If /experiment directory already exists exit and return error message and
    the full path of the present directory to the user"""
    # Check that experiment name has no special characters or spaces
    experiment = convert_to_filename(settings.experiment)
    output_path = settings.path

    # ERROR MESSAGES
    error_directory_exists = (
        f"Directory {experiment} already exists. "
        f"Snakemake will create files that are missing."
    )

    try:
        output_path.mkdir(parents=True, exist_ok=False)
        logger.info(f"Make directory: results/{experiment}")

        if settings.process_reads:
            # Raw data dir
            output_path.joinpath("raw_data").mkdir()
            logger.info(f"Make directory: results/{experiment}/raw_data/")
        if settings.parse_genebank:
            # Genome dir
            output_path.joinpath("genome_lookup").mkdir()
            logger.info(f"Make directory: results/{experiment}/genome_lookup/")
    except FileExistsError:
        logger.info(error_directory_exists)
        return


def convert_to_filename(sample_name):
    """
    Convert to a valid filename.

    Removes leading/trailing whitespace, converts internal spaces to underscores.
    Allows only alphanumeric, dashes, underscores, unicode.
    """
    return re.sub(r"(?u)[^-\w]", "", sample_name.strip().replace(" ", "_"))


def count_lines(filename: str, fastq=False):
    """
    Returns number of lines in given file for progress bar.
    Uses subprocess to call wc in shell
    """
    try:
        # Check if gz compressed
        if filename.split(".")[-1] in ["gz"]:
            # Build pipe for uncompressed file
            pipe = sub.Popen(("gzcat", "-cd", filename), stdout=sub.PIPE)
            output = int(sub.check_output(["wc", "-l"], stdin=pipe.stdout))
        # Other files...
        else:
            output = int(
                sub.check_output(["wc", "-l", filename], stderr=sub.DEVNULL)
                .split()[0]
                .decode("UTF-8")
            )
        if fastq:
            return int(output / 4)
        else:
            return int(output)
    except sub.CalledProcessError as e:
        # Catching error returns arbitrary value
        logger.debug(f"{str(e)}: using arbitrary limit for progress bar")
        return None


def read_gene_file(gene_file):
    """
    Helper for reading a sample gene file.
    Returns a dictionary where key is gene and value is cpm

    :param gene_file:
    :return:
    """
    gene_dict = dict()
    with open(gene_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for record in reader:
            gene_dict[record["locus_tag"]] = [float(record["cpm"])]

    return gene_dict


def get_version():
    """
    Returns pyinseq version
    :return:
    """
    dir = os.path.dirname(__file__)
    with open(os.path.join(dir, "VERSION"), "r") as f:
        return f.read()


def get_config_dict(config_file):
    """Loads data from config_file into dict"""
    with open(config_file, "r") as f:
        return yaml.load(f, Loader=yaml.FullLoader)


def write_config_file(args=None, format="yml", default=False) -> Path:
    """

    Given pyinseq arguments, creates a config file that can be given to snakemake

    :param args: Namespace object with pyinseq args
    :param format: file extension for config_file (json or yaml)
    :param default: boolean to write default config file
    :return: path to config_file
    """
    args_dict = vars(args)
    dumper = {"yml": yaml.dump, "json": json.dump}[format]
    if default:
        logger.info("Writing DEFAULT configuration file as template for pyinseq")
        with open(f"default-config-pyinseq.{format}", "w") as f:
            config_dict = {key: args_dict[key] for key in args_dict}
            dumper(config_dict, f)
        logger.info("Exiting gracefully...")
        exit(0)
    else:
        logger.info("Writing configuration file from provided arguments")
        config_file = Path(f"{args.experiment}-config.{format}")
        with open(config_file, "w") as f:
            if not args_dict["command"]:
                args_dict["command"] = "pyinseq"
            config_dict = {key: args_dict[key] for key in args_dict}
            dumper(config_dict, f)
        return config_file


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


"""

Functions that handle terminal operations

"""


def copy_file(src, dest):
    """Literally does what it says..."""
    shutil.copy(src, dest)


def execute(cmd) -> None:
    """Execute snakemake workflow using subprocess module"""
    # Call bash command, suppress pyinseq output
    if not isinstance(cmd, list):
        cmd = cmd.split(" ")
    process = sub.Popen(cmd)
    try:
        process.wait()
    except KeyboardInterrupt:
        # Kill if CTRL+C event occurs
        process.kill()
    finally:
        process.terminate()


def pickle_object(object, pickle_output):
    """Serializes objects that can be picked up by other processes"""
    with open(pickle_output, "wb") as f:
        pickle.dump(object, f)


def load_pickle(byte_file):
    """Loads a pickled object"""
    with open(byte_file, "rb") as f:
        return pickle.load(f)
