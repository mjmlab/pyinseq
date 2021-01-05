#!/usr/bin/env python3

import io
import re
import sys
import tqdm
import logging
import subprocess
from pathlib import Path


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
            output_path.joinpath('raw_data').mkdir()
            logger.info(f"Make directory: results/{experiment}/raw_data/fdefrefer")
        if settings.parse_genbank_file:
            # Genome dir
            output_path.joinpath('genome_lookup').mkdir()
            logger.info(f"Make directory: results/{experiment}/genome_lookup/")
    except FileExistsError:
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
        logger.debug(f"{str(e)}: using arbitrary limit for progress bar")
        return 10e10


# ===== Start here ===== #
def main():
    pass


if __name__ == "__main__":
    main()
