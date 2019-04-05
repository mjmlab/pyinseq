#!/usr/bin/env python3

import logging
import os
import re
import sys
import io
import subprocess
import tqdm
import gzip


class TqdmStream(io.StringIO):
    """Custom stream with logging compatibility for tqdm progress bars."""

    def __init__(self):
        # Call parent constructor for instance variables
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


def count_sequences(filename):
    """
    Returns number of sequences in given file for progress bar.

    Uses subprocess calls for efficiency
    """
    try:
        # Check if gz compressed
        if filename.split(".")[-1] in ["gz"]:
            # Build pipe for uncompressed file
            pipe = subprocess.Popen(("gzcat", "-d", filename), stdout=subprocess.PIPE)
            output = int(subprocess.check_output(["wc", "-l"], stdin=pipe.stdout))
        # If fastq format
        elif filename.split(".")[-1] in ["fastq"]:
            output = int(
                subprocess.check_output(
                    ["wc", "-l", filename], stderr=subprocess.DEVNULL
                )
                .split()[0]
                .decode("UTF-8")
            )
        return int(output / 4)
    except subprocess.CalledProcessError as e:
        # Catching error returns arbitrary value
        logger.warning(f"{str(e)}: using arbitrary limit for progress bar")
        return 10e10


# ===== Start here ===== #


def main():
    pass


if __name__ == "__main__":
    main()
