#!/usr/bin/env python
"""
Build bowtie index for specified genome

Specify exact format of input .fna and .ptt files here

What about RNA genes?

"""

import subprocess
import config # config file

def bowtie_build():
    organism = 'k'
    fna = organism + '.fna'
    print('\n===== Building bowtie index files for organism {} =====\n'.format(organism))
    subprocess.check_call([config.bowtie_build, '-q', fna, organism])

def bowtie_map():
    organism = 'k'
    fna = organism + '.fna'
    reads = 'Exp001_assigned.fasta'
    bowtie_output = 'bowtie_output.txt'
    print('\n===== Mapping reads to bowtie index for organism {} =====\n'.format(organism))
    subprocess.check_call([config.bowtie, '-m 1 --best --strata -a --fullref -n 1 -l 17',
        organism, '-f', reads, bowtie_output])


# ===== Start here ===== #

def main():
    bowtie_map()

if __name__ == "__main__":
    main()
