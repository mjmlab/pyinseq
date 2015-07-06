#!/usr/bin/env python
"""
Build bowtie index for specified genome

Run bowtie

(Will this work in a defined subdirectory - name it based on Experiment???)

"""

import subprocess
import config # config file
import sys # sys.argV

def bowtie_build():
    organism = sys.argv[1]
    fna = organism + '.fna'
    print('\n===== Building bowtie index files for organism {} =====\n'.format(organism))
    subprocess.check_call([config.bowtie_build, '-q', fna, organism])

def bowtie_map():
    organism = sys.argv[1]
    fna = organism + '.fna'
    reads = sys.argv[2]
    bowtie_output = 'bowtie_output.txt'
    print('\n===== Mapping reads to bowtie index for organism {} =====\n'.format(organism))
    subprocess.check_call([config.bowtie, '-m 1 --best --strata -a --fullref -n 1 -l 17',
        organism, '-q', reads, bowtie_output])


# ===== Start here ===== #

def main():
    bowtie_map()

if __name__ == "__main__":
    main()
