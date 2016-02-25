"""Configuration file for pyinseq.py package."""
import os
import sys

# TRANSPOSON END SEQUENCES
# Typical mariner transposon ends for INSeq are the defaults (same left/right ends)
# This sequence is filtered from the reads prior to mapping
# Mariner transposons inserts at TA dinucleotides; do not include the corresponding
# TA from the tranpsoson here
transposonLeft = 'ACAGGTTG'
transposonRight = 'ACAGGTTG'

# Bowtie package selected (linux/mac/windows) based on operating system detected
platform = sys.platform.lower()
current_folder = os.path.dirname(os.path.abspath(__file__))
packages_folder = os.path.join(current_folder, '..', 'packages')

if sys.platform.lower().startswith('linux'):
    # linux bowtie path
    bowtie = os.path.join(packages_folder, 'bowtie-1.1.1-linux', 'bowtie')

elif sys.platform.lower().startswith('darwin'):
    # mac bowtie path
    bowtie = os.path.join(packages_folder, 'bowtie-1.1.1-mac', 'bowtie')

elif sys.platform.lower().startswith('win'):
    # windows bowtie path
    bowtie = os.path.join(packages_folder, 'bowtie-1.1.1-win', 'bowtie')

# PATH TO BOWTIE-BUILD (appends '-build' on the path above)
bowtieBuild = '{bowtiepath}-build'.format(bowtiepath=bowtie)
