#!/usr/bin/env python
"""
Build bowtie index for specified genome

Run bowtie

"""

import subprocess
import config # config file

# organism is the name of the .fna file (can be multifasta)
# organism also becomes the base name for the bowtie index
def bowtie_build(organism):
    fna = organism + '.fna'
    print('Building bowtie index files for organism {}'.format(organism))
    subprocess.check_call([config.bowtie_build, '-q', fna, organism])

# organism is the base name for the bowtie index
def bowtie_map(organism, reads, bowtie_output):
    fna = organism + '.fna'
    #bowtie_output = 'bowtie_output.txt'
    print('Mapping reads to bowtie index for organism {}'.format(organism))
    # String version of the shell command
    bash_command = '{0} -m 1 --best --strata -a --fullref -n 1 -l 17 {1} -q {2} {3} -p 2' \
        .format(config.bowtie, organism, reads, bowtie_output)
    # Convert bash command to run properly - no spaces; instead list of entries
    # that will appear in the shell as space-separated
    subprocess.check_call(bash_command.split(' '))

# ===== Start here ===== #

def main():
    pass

if __name__ == '__main__':
    main()
