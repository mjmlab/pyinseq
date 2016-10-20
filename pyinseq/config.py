#!/usr/bin/env python3

'''Configuration file for pyinseq.py package.'''
import os
import sys

# TRANSPOSON END SEQUENCES
# Typical mariner transposon ends for INSeq are the defaults (same left/right ends)
# This sequence is filtered from the reads prior to mapping
# Mariner transposons inserts at TA dinucleotides; do not include the corresponding
# TA from the tranpsoson here
transposon_end = {'left': 'ACAGGTTG', 'right': 'ACAGGTTG'}

# Bowtie package selected (linux/mac/windows) based on operating system detected
platform = sys.platform.lower()
current_folder = os.path.dirname(os.path.abspath(__file__))
packages_folder = os.path.join(current_folder, 'third-party')

default_version = 'bowtie-1.1.1-mac'

if platform.startswith('linux'):
    # linux bowtie path
    default_version = 'bowtie-1.1.1-linux'

elif platform.startswith('win'):
    # windows bowtie path
    default_version = 'bowtie-1.1.1-win'

bowtie = os.path.join(packages_folder, default_version, 'bowtie')

# PATH TO BOWTIE-BUILD (appends '-build' on the path above)
bowtieBuild = '{bowtiepath}-build'.format(bowtiepath=bowtie)
