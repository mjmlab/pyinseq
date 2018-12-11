#!/usr/bin/env python3

import logging
import re
import subprocess
import pyinseq.config as config

logger = logging.getLogger("pyinseq")


def bowtie_build(organism):
    """Build a bowtie index given a fasta nucleotide file."""
    fna = organism + ".fna"
    subprocess.check_call([config.bowtieBuild, "-q", fna, organism])


def bowtie_map(organism, reads, bowtie_output):
    """Map fastq reads to a bowtie index."""
    fna = organism + ".fna"
    # String version of the shell command
    bash_command = f"{config.bowtie} -m 1 --best --strata -a --fullref -n 1 -l 17 {organism} -q {reads} {bowtie_output} -p 2"
    # Convert bash command to run properly - no spaces; instead list of entries
    # that will appear in the shell as space-separated
    # Consider shlex.split() instead of split() -- any benefit here?
    # subprocess.check_call(bashCommand.split(' '))
    proc = subprocess.Popen(
        bash_command.split(" "), stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    bowtie_msg_out = proc.communicate()[0]
    # decode from bytes to unicode
    bowtie_msg_out = bowtie_msg_out.decode("utf8")
    return bowtie_msg_out


def parse_bowtie(bowtie_message):
    """Parse bowtie results into a dictionary."""
    bowtie_msg_dict = {}
    for line in bowtie_message.split("\n"):
        # extract counts from bowtie printing
        m = re.search("^(\#.+:) (\d+)", line)
        try:
            message, count = m.group(1), int(m.group(2))
            bowtie_msg_dict[message] = count
        except:
            pass
    return bowtie_msg_dict


# ===== Start here ===== #


def main():
    pass


if __name__ == "__main__":
    main()
