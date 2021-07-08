#!/usr/bin/env python3

"""

Contains functions for parsing arguments from command line

"""

import os
import argparse

# Module imports
from pyinseq.utils import get_version


def get_snake_parser():
    # Snake parser will be parent for all other parsers
    snake_parser = argparse.ArgumentParser(add_help=False)
    snake_group = snake_parser.add_argument_group("SNAKEMAKE")
    snake_group.add_argument(
        "--get_default_config",
        action="store_true",
        help="Writes a default configuration file for pyinseq run that can be modified",
    )
    snake_group.add_argument(
        "--config_format",
        default="yml",
        help="Format for the configuration file. Options are: 'json' or 'yml'",
    )
    snake_group.add_argument(
        "-c",
        "--config",
        help="Provide config file (in json or yaml format) instead of pyinseq arguments to run pipeline",
        required=False,
    )
    snake_group.add_argument(
        "-t",
        "--threads",
        default=os.cpu_count(),
        help="Number of threads that snakemake can use to run pyinseq, the more the better for parallel processing",
        type=int,
    )
    snake_group.add_argument(
        "--snakemake_params",
        help="Additional params passed to snakemake. "
        "This should be included at the end of the command since everything will be passed to snakemake so make sure they are correct. "
        "For example, you can use `-n` to check which files snakemake will create without execution the full workflow.",
        nargs="...",
        default=[],
        type=str,
    )
    snake_group.add_argument(
        "--test",
        action="store_true",
        help="Runs pytest on pyinseq module",
        required=False,
        default=False,
    )
    return snake_parser


def get_args():
    """Parse command line arguments for main pyinseq. Returns both parser and Namespace object"""

    parser = argparse.ArgumentParser("pyinseq", parents=[get_snake_parser()])
    parser.add_argument(
        "-v", "--version", action="version", version=f"pyinseq: {get_version()}"
    )
    parser.add_argument(
        "-i", "--input", help="input Illumina reads file or folder", required=False
    )
    parser.add_argument(
        "-s", "--samples", help="sample list with barcodes", required=False
    )
    parser.add_argument(
        "-e",
        "--experiment",
        help="experiment name (no spaces or special characters)",
        required=False,
    )
    parser.add_argument(
        "-g",
        "--genome",
        help="genome in GenBank format (one concatenated file for multiple contigs/chromosomes)",
        required=False,
    )
    parser.add_argument(
        "-d",
        "--disruption",
        help="fraction of gene disrupted (0.0 - 1.0). Default: 0.9",
        default=0.9,
        type=float,
    )
    parser.add_argument(
        "--min_count",
        help="Minimum number of reads per insertion site",
        default=3,
        type=int,
    )
    parser.add_argument(
        "--max_ratio",
        help="Maximum ratio of left:right or right:left reads per insertion site",
        default=10,
        type=int,
    )
    parser.add_argument(
        "--barcode_length",
        help="Length of the barcode which is used to demultiplex samples (4 - 16bp)",
        type=int,
        default=4,
    )
    parser.add_argument(
        "--transposon_seq",
        help="Sequence for the transposon that flanks reads",
        default="ACAGGTTG",
    )
    parser.add_argument(
        "--gff3",
        help="generate GFF3 file",
        action="store_true",
        required=False,
        default=False,
    )
    subparsers = parser.add_subparsers(
        dest="command", help="Command for performing specialized tasks"
    )
    # demultiplex
    sub_parser_demultiplex = subparsers.add_parser(
        "demultiplex",
        parents=[get_snake_parser()],
        help="Demultiplex reads into barcode samples",
    )
    sub_parser_demultiplex.add_argument(
        "-i", "--input", help="input Illumina reads file or folder", required=False
    )
    sub_parser_demultiplex.add_argument(
        "-s", "--samples", help="sample list with barcodes", required=False
    )
    sub_parser_demultiplex.add_argument(
        "-e",
        "--experiment",
        help="experiment name (no spaces or special characters)",
        required=False,
    )
    sub_parser_demultiplex.add_argument(
        "--notrim",
        help="do not write trimmed reads (i.e. write raw reads only)",
        action="store_true",
        required=False,
        default=False,
    )
    sub_parser_demultiplex.add_argument(
        "--barcode_length",
        help="Length of the barcode which is used to demultiplex samples (4 - 16bp)",
        type=int,
        default=4,
    )
    sub_parser_demultiplex.add_argument(
        "--transposon_seq",
        help="Sequence for the transposon that flanks reads",
        default="ACAGGTTG",
    )

    # Genomeprep
    sub_parser_genome_prep = subparsers.add_parser(
        "genomeprep",
        parents=[get_snake_parser()],
        help="Prepare genome files from nucleotide sequences",
    )
    sub_parser_genome_prep.add_argument(
        "-e",
        "--experiment",
        help="experiment name (no spaces or special characters)",
        required=False,
    )
    sub_parser_genome_prep.add_argument(
        "-g",
        "--genome",
        help="genome in GenBank format (one concatenated file for multiple contigs/chromosomes)",
        required=False,
    )
    sub_parser_genome_prep.add_argument(
        "--noindex",
        help="do not generate bowtie indexes",
        action="store_true",
        required=False,
        default=False,
    )
    sub_parser_genome_prep.add_argument(
        "--gff3",
        help="generate GFF3 file",
        action="store_true",
        required=False,
        default=False,
    )
    return parser, parser.parse_args()


"""Inactive arguments in current version
parser.add_argument('--nobarcodes',
                    help='barcodes have already been removed from the samples; \
                    -i should list the directory with filenames (.fastq.gz) \
                    corresponding to the sample names',
                    action='store_true',
                    default=False)
parser.add_argument('--compress',
                    help='compress (gzip) demultiplexed samples',
                    action='store_true',
                    default=False)
parser.add_argument('--keep_all',
                    help='keep all intermediate files generated',
                    action='store_true',
                    default=False)"""
