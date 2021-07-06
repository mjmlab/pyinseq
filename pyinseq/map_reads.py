#!/usr/bin/env python3

"""

Functions for using bowtie

"""

import re
import logging
import subprocess

# Module imports
from pyinseq.logger import pyinseq_logger, add_fileHandler

# Custom snake message for jobs
SNAKE_MESSAGE = """
(SNAKEMAKE INFO)
Job counts:
\tcount	jobs
\t1\tbowtie_mapping
\t1
"""

formatter = logging.Formatter(
    fmt="{asctime} - {levelname} - map_reads - {message}",
    datefmt="%Y-%m-%d %H:%M",
    style="{",
)


def bowtie_build(organism):
    """Build a bowtie index given a fasta nucleotide file."""
    fna = organism + ".fna"
    subprocess.check_call(["bowtie-build", "-q", fna, organism])


def bowtie_map(organism, reads, bowtie_output, threads):
    """Map fastq reads to a bowtie index."""
    fna = organism + ".fna"
    # String version of the shell command
    bash_command = (
        f"bowtie -m 1 --best --strata -a --fullref -n 1 -l 17 "
        f"{organism} -q {reads} {bowtie_output} -p {threads}"
    )
    # Convert bash command to run properly - no spaces; instead list of entries
    # that will appear in the shell as space-separated
    # Consider shlex.split() instead of split() -- any benefit here?
    # subprocess.check_call(bashCommand.split(' '))
    logger.info(f"Mapping reads in {reads} to genome using bowtie")
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
        m = re.search(r"^(\#.+:) (\d+)", line)
        try:
            message, count = m.group(1), int(m.group(2))
            bowtie_msg_dict[message] = count
        except:
            pass
    return bowtie_msg_dict


def summarize_mapping(sample, bowtie_msg_dict):
    """Summarizes bowtie mapping step into pyinseq_logger io"""
    pyinseq_logger.logger_io.write(
        f"- Mapped reads from {sample} to genome\nBOWTIE OUTPUT - {sample}:\n"
    )
    for key, value in bowtie_msg_dict.items():
        pyinseq_logger.logger_io.write(f"{key} {value}")
    return


# ===== Start here ===== #
if __name__ == "__main__":
    # Modified to work with snakemake
    organism = snakemake.params[0]
    reads = snakemake.input[0]
    output = snakemake.output[0]
    threads = snakemake.threads
    sample = snakemake.wildcards[0]

    # Setup logger
    pyinseq_logger.setup_logger(formatter=formatter)
    logger = pyinseq_logger.logger
    pyinseq_logger.summary_log = snakemake.params.summary_log
    add_fileHandler(logger, snakemake.params.log, formatter=formatter)

    # Print and write custom snake message
    print(SNAKE_MESSAGE)
    pyinseq_logger.snake_io.write(SNAKE_MESSAGE)

    # Begin bowtie
    bowtie_msg = bowtie_map(organism, reads, output, threads)
    # Print bowtie results for sample
    logger.info(f"Bowtie results for {sample}:")
    with open(snakemake.params.log, "a+") as f:
        print(bowtie_msg)
        print(bowtie_msg, file=f)

    bowtie_msg_dict = parse_bowtie(bowtie_msg)
    summarize_mapping(sample, bowtie_msg_dict)
    # Writes to log files
    pyinseq_logger.summarize_step()
