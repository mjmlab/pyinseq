#!/usr/bin/env python

import os

def createSamplesDirectory():
    """
    Create the samples directory if needed
    """

    # Create /pyinseq/samples/ path if does not exist already
    try:
        os.makedirs('samples/')
    except:
        pass

def createExperimentDirectories(experiment):
    """
    Create the project directory and subdirectories

    Attempt to create the directory structure:
    /[experiment]
        /temp

    If /experiment directory already exists exit and return error message and
    the full path of the present directory to the user"""

    # Check that experiment name has no special characters or spaces
    pass

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
