#!/usr/bin/env python3
"""
Demultiplexes a FASTQ file into multiple files by 5' end barcode.

Output path includes the experiment and sample name:
(pyinseq)/results/{experiment}/{sample}.fastq

"""
import tqdm
import logging
import re
import screed
import sys
import os
from .utils import count_sequences

logger = logging.getLogger("pyinseq")


def demultiplex_fastq(reads, samples_dict: dict, settings) -> int:
    """Demultiplex a fastq input file by 5' barcode into separate files.

    Use regex to identify the chromosome slice and save this in the read
    record as 'trim', e.g., read['trim'] = (4, 21)
    Save raw reads into '{experiment}/raw_data/{sample_name}.fastq'
    Save trimmed reads into '{experiment}/{sample_name}_trimmed.fastq'
    """
    # Dictionary of lists to hold FASTQ reads until they are written to files
    # Keys are barcodes
    demultiplex_dict = {}
    for sample in samples_dict:
        bc = samples_dict[sample]["barcode"]
        demultiplex_dict[bc] = []
    demultiplex_dict["other"] = []  # unassigned barcodes
    # count of reads
    n_reads = 0
    # For each line in the FASTQ file:
    #   Assign to barcode (fastq record into the dictionary)
    #   Cache by barcode; write to the appropriate output files (untrimmed and trimmed)
    #   Use triple brackets in f-strings to leave a single bracket for regex
    # REGEX DESCRIPTION:
    #   beginning of string
    #   group(1) = barcode, any 4-bp of mixed ACGT
    #   group(2) = 16-17 bp of chromosmal sequence
    #   flanking transposon sequence, last two must be TA for transposon

    pattern = re.compile(
        f"""
        ^                                     #   beginning of string
        ([ACGT]{{{settings.barcode_length}}}) #   group(1) = barcode, any bp (4-16bp) of mixed ACGT
        ([NACGT][ACGT]{{13,14}}               #   group(2) = 16-17 bp of chromosmal sequence
        (?:TA)){settings.transposon_seq}""",  #   flanking transposon sequence, last two must be TA for transposon
        re.VERBOSE,

    )

    # Initiate progress bar
    logger.info("Preparing to demultiplex reads")

    p_bar = tqdm.tqdm(
        total=count_sequences(reads),
        desc=f"Demultiplexing Reads",
        unit="reads",
        leave=False,
    )
    for read in screed.open(reads):
        p_bar.update()
        m = re.search(pattern, read.sequence)
        try:
            barcode, chrom_seq = m.group(1), m.group(2)
            if barcode in demultiplex_dict:
                read["trim"] = m.span(2)  # trim slice
                demultiplex_dict[barcode].append(read)
            else:
                demultiplex_dict["other"].append(read)
        except:
            demultiplex_dict["other"].append(read)
        # Every 10^6 sequences write and clear the dictionary
        n_reads += 1
        if n_reads % 5e6 == 0:
            # Disable Progress Bar
            p_bar.close()
            logger.info(
                "Demultiplexed {:,} samples and writing trimmed reads".format(n_reads)
            )
            # Enable Progress Bar
            p_bar.disable = False
            p_bar.refresh()
            write_reads(demultiplex_dict, samples_dict, settings)
            # Write trimmed reads only when needed
            if settings.write_trimmed_reads:
                write_trimmed_reads(demultiplex_dict, samples_dict, settings)
            # Clear the dictionary after writing to file
            for sample_name in demultiplex_dict:
                demultiplex_dict[sample_name] = []
    # Remove any traces of Progress Bar
    p_bar.close()
    write_reads(demultiplex_dict, samples_dict, settings)
    # Write trimmed reads only when needed
    if settings.write_trimmed_reads:
        write_trimmed_reads(demultiplex_dict, samples_dict, settings)

    logger.info(f"Total records demultiplexed: {n_reads:,}")
    return n_reads


def write_reads(demultiplex_dict, samples_dict, settings):
    """Write the fastq data to the correct (demultiplexed) file."""
    # inverting the value: key pairs to barcode: sample, and also adding 'other': '_other'
    barcode_dict = {"other": "_other"}
    for sample in samples_dict:
        barcode_dict[samples_dict[sample]["barcode"]] = sample
    for barcode in demultiplex_dict:
        if demultiplex_dict[barcode]:
            with open(
                f"{settings.path}raw_data/{barcode_dict[barcode]}.fastq", "a"
            ) as fo:
                logger.debug(f"Writing reads to {barcode_dict[barcode]}.fastq")
                # Initialize Progress Bar
                p_bar = tqdm.tqdm(
                    demultiplex_dict[barcode],
                    total=len(demultiplex_dict[barcode]),
                    desc=f"Writing {barcode_dict[barcode]} reads",
                    unit="reads",
                    leave=False,
                    position=1,
                )
                for read in demultiplex_dict[barcode]:
                    p_bar.update()
                    fo.write(f"@{read.name}\n{read.sequence}\n+\n{read.quality}\n")
                # Remove any traces of Progress Bar
                p_bar.close()
    return


def write_trimmed_reads(demultiplex_dict, samples_dict, settings):
    """Write the fastq data to the correct (demultiplexed) file."""
    # inverting the value: key pairs to barcode: sample. Exclude 'other' here
    barcode_dict = {}
    for sample in samples_dict:
        barcode_dict[samples_dict[sample]["barcode"]] = sample
    for barcode in demultiplex_dict:
        if barcode != "other":
            with open(
                f"{settings.path}/{barcode_dict[barcode]}_trimmed.fastq", "a"
            ) as fo:
                logger.debug(f"Writing trimmed reads to {barcode_dict[barcode]}.fastq")
                # Initialize Progress Bar
                p_bar = tqdm.tqdm(
                    demultiplex_dict[barcode],
                    total=len(demultiplex_dict[barcode]),
                    desc=f"Writing {barcode_dict[barcode]} trimmed reads",
                    unit="reads",
                    leave=False,
                    position=1,
                )
                for read in demultiplex_dict[barcode]:
                    p_bar.update()
                    fo.write(
                        f"@{read.name}\n"
                        f"{read.sequence[slice(read.trim[0], read.trim[1])]}\n"
                        f"+\n"
                        f"{read.quality[slice(read.trim[0], read.trim[1])]}\n"
                    )
                # Remove any traces of Progress Bar
                p_bar.close()
    return


def main():
    """Start here."""


if __name__ == "__main__":
    main()
