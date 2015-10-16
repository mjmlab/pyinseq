#!/usr/bin/env python

import os
import re

def createSamplesDirectory():
    """
    Create the samples directory if needed
    """

    # Create /pyinseq/samples/ path if does not exist already
    try:
        os.makedirs('samples/')
    except:
        pass

def createDemultiplexFiles(experiment, sample_name):
    """
    Create blank samples/{experiment}/{sample_name}.fastq files
    """

    #TODO: Future, ask if user wants files overwritten?

    # Create /pyinseq/samples/experiment/sample_name.fastq if does not exist already
    directory = 'samples/{0}/'.format(experiment)
    filename = '{0}.fastq.gz'.format(sample_name)
    filepath = '{0}{1}'.format(directory, filename)

    try:
        os.makedirs(directory)
    except OSError:
        pass

    if not os.path.exists(filepath):
        open(filepath, 'w').close()
    else:
        print('Error: File {} exists already.'.format(filepath))
        exit(1)


def createExperimentDirectories(experiment):
    """
    Create the project directory and subdirectories

    Attempt to create the directory structure:
    /[experiment]
        /temp

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
        os.makedirs('{}/temp/'.format(experiment))
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
    #args = parseArgs(sys.argv[1:])
    #print(args.input)
    #print(args.samples)
    #assignAndTrim(args.input, args.samples)
    experiment = 'E1'

    createDirectories(experiment)

if __name__ == '__main__':
    main()
