#!/usr/bin/env python
"""
Build bowtie index for specified genome

Run bowtie

"""

import config # config file
import subprocess

# organism is the name of the .fna file (can be multifasta)
# organism also becomes the base name for the bowtie index
def bowtieBuild(organism):
    fna = organism + '.fna'
    subprocess.check_call([config.bowtieBuild, '-q', fna, organism])

# organism is the base name for the bowtie index
def bowtieMap(organism, reads, bowtieOutput):
    fna = organism + '.fna'
    # String version of the shell command
    bashCommand = '{0} -m 1 --best --strata -a --fullref -n 1 -l 17 {1} -q {2} {3} -p 2' \
        .format(config.bowtie, organism, reads, bowtieOutput)
    # Convert bash command to run properly - no spaces; instead list of entries
    # that will appear in the shell as space-separated
    # Consider shlex.split() instead of split() -- any benefit here?
    subprocess.check_call(bashCommand.split(' '))

# ===== Start here ===== #

def main():
    pass

if __name__ == '__main__':
    main()
