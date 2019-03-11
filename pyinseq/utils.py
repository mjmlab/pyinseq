#!/usr/bin/env python3

import logging
import mmap
import os
import re
import sys
import io
import screed.fastq
import tqdm


class TqdmStream(io.StringIO):
    """Custom stream for logging compatibility with tqdm progress bars."""

    def __init__(self):
        io.StringIO.__init__(self)

    def write(self, s):
        # Redirect through tqdm.write and output to stderr
        tqdm.tqdm.write(s, file=sys.stderr, end="")


# This controls the stdout logging.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(module)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M",
    stream=TqdmStream(),
)

logger = logging.getLogger("pyinseq")


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

    # ERROR MESSAGES
    error_directory_exists = (
        "PyINSeq Error: The directory already exists for experiment {0}\n"
        f"Delete or rename the {experiment} directory, or provide a new experiment\n"
        "name for the current analysis"
    )

    # Create path or exit with error if it exists.
    try:
        if settings.process_reads:
            os.makedirs(f"results/{experiment}/raw_data/")
            logger.info(f"Make directory: results/{experiment}")
            logger.info(f"Make directory: results/{experiment}/raw_data/")
        # Only make the genome lookup directory if needed
        if settings.parse_genbank_file:
            os.makedirs(f"results/{experiment}/genome_lookup/")
            logger.info(f"Make directory: results/{experiment}/genome_lookup/")

    # TODO: Make analysis directory using path attributes from Settings
    except OSError:
        print(error_directory_exists)
        exit(1)


def convert_to_filename(sample_name):
    """
    Convert to a valid filename.

    Removes leading/trailing whitespace, converts internal spaces to underscores.
    Allows only alphanumeric, dashes, underscores, unicode.
    """
    return re.sub(r"(?u)[^-\w]", "", sample_name.strip().replace(" ", "_"))


def fastq_generator(filename):
    """ 
    Returns number of reads in given fastq file.

    Also returns a generator that yields Screed Records which can
    be looped through.

    """
    # Open file and map to mmap
    f = open(filename, "r")
    fastq_mmap = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
    # Get iterable for counting reads
    record_iterable = screed.fastq.fastq_iter(fastq_mmap)
    # Count number of Screed Records
    reads = 0
    for i in screed.fastq.fastq_iter(mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)):
        reads += 1
    return {
        "total_reads": reads,
        "reads_generator": screed.fastq.fastq_iter(fastq_mmap),
    }


# ===== Start here ===== #


def main():
    pass


if __name__ == "__main__":
    main()
