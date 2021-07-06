#!/usr/bin/env python3

"""
Demultiplexes a FASTQ file into multiple files by 5' end barcode.

Output path includes the experiment and sample name:
(pyinseq)/results/{experiment}/{sample}.fastq

"""

import re
import screed

# Module imports
from pyinseq.utils import count_lines
from pyinseq.logger import pyinseq_logger, TqdmBarLogger

logger = pyinseq_logger.logger


def demultiplex_fastq(reads, samples_dict: dict, settings) -> int:
    """Demultiplex a fastq input file by 5' barcode into separate files.

    Use regex to identify the chromosome slice and save this in the read
    record as 'trim', e.g., read['trim'] = (4, 21)
    Save raw reads into '{experiment}/raw_data/{sample_name}.fastq'
    Save trimmed reads into '{experiment}/{sample_name}_trimmed.fastq'

    :reads: path to file with reads
    :samples_dict: dict where sample is key and value={"barcode": "NNNN"}
    :settings: Settings object for experiment
    :return: a dict where key=barcode and value=num of reads
    """
    # Dictionary of lists to hold FASTQ reads until they are written to files
    # Keys are barcodes
    demultiplex_dict = {}
    # Keys are barcodes but used to count reads in sample
    demultiplex_count_dict = {}
    for sample in samples_dict:
        bc = samples_dict[sample]["barcode"]
        demultiplex_dict[bc] = []
        demultiplex_count_dict[bc] = 0
    demultiplex_dict["other"] = []  # unassigned barcodes
    demultiplex_count_dict["other"] = 0
    # For each line in the FASTQ file:
    #   Assign to barcode (fastq record into the dictionary)
    #   Cache by barcode; write to the appropriate output files (untrimmed and trimmed)
    #   Use triple brackets in f-strings to leave a single bracket for regex
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
    num_fastq_lines = count_lines(reads, fastq=True)
    tqdm_logger = TqdmBarLogger(
        logger, num_fastq_lines, desc="Demultiplexing reads", unit="reads"
    )
    # count of reads
    n_reads = 0
    for read in screed.open(reads):
        tqdm_logger.update()
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
        # Buffer dumping reads into barcode file
        if n_reads % 5e6 == 0:
            tqdm_logger.info(f"Demultiplexed {n_reads:,} reads. Writing to file")
            write_reads(demultiplex_dict, samples_dict, settings)
            # Clear the dictionary after writing to file
            for barcode in demultiplex_dict:
                # Record amount of reads in barcode
                demultiplex_count_dict[barcode] += len(demultiplex_dict[barcode])
                demultiplex_dict[barcode] = []
    write_reads(demultiplex_dict, samples_dict, settings)
    # Clear the dictionary after writing to file
    for barcode in demultiplex_dict:
        # Record amount of reads in barcode
        demultiplex_count_dict[barcode] += len(demultiplex_dict[barcode])
        demultiplex_dict[barcode] = []
    # Remove any traces of Progress Bar
    tqdm_logger.p_bar.close()
    # print demultiplex summary
    pyinseq_logger.logger_io.write(f"{n_reads} reads demultiplexed\n")
    for bc in demultiplex_count_dict:
        for sample in samples_dict:
            if samples_dict[sample]["barcode"] == bc:
                samples_dict[sample]["demultiplexed_reads"] = demultiplex_count_dict[bc]
                pyinseq_logger.logger_io.write(
                    f"- {sample} ({bc}): {demultiplex_count_dict[bc]:,} reads"
                )
    logger.info(f"Total reads demultiplexed: {n_reads:,}")

    # Return an update samples_dict
    return samples_dict


def write_reads(demultiplex_dict, samples_dict, settings):
    """Write the fastq data to the correct (demultiplexed) file."""
    # inverting the value: key pairs to barcode: sample, and also adding 'other': '_other'
    barcode_dict = {"other": "_other"}
    for sample in samples_dict:
        barcode_dict[samples_dict[sample]["barcode"]] = sample
    # Write to sample file
    for barcode in demultiplex_dict:
        if demultiplex_dict[barcode]:
            with open(
                settings.raw_path.joinpath(f"{barcode_dict[barcode]}.fastq"), "a"
            ) as fo:
                logger.debug(f"Writing reads to {barcode_dict[barcode]}.fastq")
                # Initialize Progress Bar
                num_lines_to_write = len(demultiplex_dict[barcode])
                tqdm_write_logger = TqdmBarLogger(
                    logger,
                    num_lines_to_write,
                    desc=f"Writing reads to {barcode_dict[barcode]}.fastq",
                    unit="reads",
                )
                for read in demultiplex_dict[barcode]:
                    tqdm_write_logger.update()
                    fo.write(f"@{read.name}\n{read.sequence}\n+\n{read.quality}\n")
                tqdm_write_logger.close()
        # Write trimmed if necessary
        if settings.write_trimmed_reads and barcode != "other":
            with open(
                settings.path.joinpath(f"{barcode_dict[barcode]}_trimmed.fastq"), "a"
            ) as fo:
                logger.debug(
                    f"Writing trimmed reads to {barcode_dict[barcode]}_trimmed.fastq"
                )
                # Same logger, TODO: could maybe wrap loop into iterator
                tqdm_write_logger = TqdmBarLogger(
                    logger,
                    num_lines_to_write,
                    desc=f"Writing trimmed reads to {barcode_dict[barcode]}_trimmed.fastq",
                    unit="reads",
                )
                for read in demultiplex_dict[barcode]:
                    tqdm_write_logger.update()
                    fo.write(
                        f"@{read.name}\n"
                        f"{read.sequence[slice(read.trim[0], read.trim[1])]}\n"
                        f"+\n"
                        f"{read.quality[slice(read.trim[0], read.trim[1])]}\n"
                    )
                tqdm_write_logger.close()
    return


def main():
    """Start here."""


if __name__ == "__main__":
    main()
