#!/usr/bin/env python

import os
import re


def createExperimentDirectories(experiment):
    """
    Create the project directory and subdirectories

    Attempt to create the directory structure:

    results/
    |
    +-{experiment}/        # User-specified experiment name
      |
      +-raw_data/          # For demultiplexed reads
      |
      +-genome_lookup/     # Genome fna and ftt files, bowtie indexes

    If /experiment directory already exists exit and return error message and
    the full path of the present directory to the user"""

    # Check that experiment name has no special characters or spaces
    experiment = convert_to_filename(experiment)

    # ERROR MESSAGES
    errorDirectoryExists = \
    'PyINSeq Error: The directory already exists for experiment {0}\n' \
    'Delete or rename the {0} directory, or provide a new experiment\n' \
    'name for the current analysis'.format(experiment)

    # Create path or exit with error if it exists.
    try:
        os.makedirs('results/{}/raw_data/'.format(experiment))
        print('Making directory: results/{}'.format(experiment))
        print('Making directory: results/{}/raw_data/'.format(experiment))
        os.makedirs('results/{}/genome_lookup/'.format(experiment))
        print('Making directory: results/{}/genome_lookup/'.format(experiment))
    except OSError:
        print(errorDirectoryExists)
        exit(1)

def convert_to_filename(sample_name):
    """
    Convert to a valid filename.

    Removes leading/trailing whitespace, converts internal spaces to underscores.
    Allows only alphanumeric, dashes, underscores, unicode.
    """
    return re.sub(r'(?u)[^-\w]', '', sample_name.strip().replace(' ', '_'))


# ===== Start here ===== #

def main():
    pass

if __name__ == '__main__':
    main()
