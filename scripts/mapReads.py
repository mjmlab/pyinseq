#!/usr/bin/env python
"""
Build bowtie index for specified genome

Run bowtie

"""

import config  # config file
import re
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
    # subprocess.check_call(bashCommand.split(' '))
    proc = subprocess.Popen(bashCommand.split(' '), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    bowtie_message = proc.communicate()[0]
    # decode from bytes to unicode
    bowtie_message = bowtie_message.decode('utf8')
    print(bowtie_message)
    return bowtie_message

def parseBowtie(bowtieMessage):
    bowtie_messages = {
        '# reads processed:': 0,
        '# reads with at least one reported alignment:': 0,
        '# reads that failed to align:': 0,
        '# reads with alignments suppressed due to -m:': 0
    }
    for line in bowtieMessage.split('\n'):
        # extract counts from bowtie printing
        m = re.search('^(\#.+:) (\d+)', line)
        try:
            message, count = m.group(1), int(m.group(2))
            bowtie_messages[message] = count
        except:
            pass


# ===== Start here ===== #

def main():
    pass

if __name__ == '__main__':
    main()
