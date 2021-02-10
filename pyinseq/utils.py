#!/usr/bin/env python3
"""

 Functions that have no home...

"""
import re
import os
import csv
import yaml
import json
import subprocess
from pathlib import Path
# Module imports
from pyinseq.logger import logger as logger


def create_experiment_directories(settings, snakemake=False):
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
    output_path = Path(f"results/{experiment}")

    # ERROR MESSAGES
    error_directory_exists = (
        f"PyINSeq Error: The directory already exists for experiment {experiment}\n"
        f"Delete or rename the {experiment} directory, or provide a new experiment\n"
        "name for the current analysis"
    )

    try:
        output_path.mkdir(parents=True, exist_ok=False)
        logger.info(f"Make directory: results/{experiment}")

        if settings.process_reads:
            # Raw data dir
            output_path.joinpath("raw_data").mkdir()
            logger.info(f"Make directory: results/{experiment}/raw_data/")
        if settings.parse_genbank_file:
            # Genome dir
            output_path.joinpath("genome_lookup").mkdir()
            logger.info(f"Make directory: results/{experiment}/genome_lookup/")
    except FileExistsError:
        if snakemake:
            logger.info(
                f"Directory {experiment} already exists. Since you decided to use snakemake, we will skip this...")
            return
        print(error_directory_exists)
        exit(1)
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
            pipe = subprocess.Popen(("gzcat", "-d", filename), stdout=subprocess.PIPE)
            output = int(subprocess.check_output(["wc", "-l"], stdin=pipe.stdout))
        # Other files...
        else:
            output = int(
                subprocess.check_output(
                    ["wc", "-l", filename], stderr=subprocess.DEVNULL
                )
                .split()[0]
                .decode("UTF-8")
            )
        if fastq:
            return int(output / 4)
        else:
            return int(output)
    except subprocess.CalledProcessError as e:
        # Catching error returns arbitrary value
        logger.debug(f"{str(e)}: using arbitrary limit for progress bar")
        return 10e10


def read_gene_file(gene_file):
    """
    Helper for reading a sample gene file.
    Returns a dictionary where key is gene and value is cpm

    :param gene_file:
    :return:
    """
    gene_dict = dict()
    with open(gene_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for record in reader:
            gene_dict[record['locus_tag']] = [float(record['cpm'])]

    return gene_dict


def get_version():
    """
    Returns pyinseq version
    :return:
    """
    dir = os.path.dirname(__file__)
    with open(os.path.join(dir, "VERSION"), 'r') as f:
        return f.read()


def get_config_dict(config_file):
    """ Loads data from config_file into dict"""
    with open(config_file, 'r') as f:
        return yaml.load(f, Loader=yaml.FullLoader)


TEMPLATE = {
    'input': None,
    'samples': None,
    'genome': None,
    'experiment': None,
    'threads': None,
}

def write_config_file(args, format='yaml', default=False) -> Path:
    """

    Given pyinseq arguments, creates a config file that can be given to snakemake

    :param args: Namespace object with pyinseq args
    :param format: file extension for configfile (json or yaml)
    :return:
    """
    dumper = {'yaml': yaml.dump, 'json': json.dump}[format]
    configfile = Path(f'{args.experiment}_pyinseq-snake_config.{format}')
    if default:
        configfile = Path(f"pyinseq-snake_default-config.{format}")
    with open(configfile, 'w') as f:
        config_dict = {key: vars(args)[key] for key in TEMPLATE}
        dumper(config_dict, f)
    return configfile

