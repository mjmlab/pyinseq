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


    """

    from Goodman:

    # mismatch replaced with 0

    $bowtie_path/bowtie -m 1 --best --strata -a --fullref -n 0 -l 17  $path\/indexes\/$index\/$index -f bcsortedseqs\/$input\_$codes_hash{$_}\.fas results\/$input\_$codes_hash{$_}\.bowtiemap\n

    """


# ===== Start here ===== #

def main():
    bowtie_build()

if __name__ == "__main__":
    main()
