#!/usr/bin/env python3

"""

Contains functions for parsing arguments from command line

"""
import os
import argparse
# Module imports
from pyinseq.utils import get_version


def get_args() -> argparse.Namespace:
    """Parse command line arguments for main pyinseq."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version', version=f"pyinseq: {get_version()}",
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
        "-d", "--disruption", help="fraction of gene disrupted (0.0 - 1.0)", default=1.0
    )
    parser.add_argument(
        "--min_count", help="Minimum number of reads per insertion site", default=3
    )
    parser.add_argument(
        "--max_ratio",
        help="Maximum ratio of left:right or right:left reads per insertion site",
        default=10,
    )
    """Inactive arguments in current version
    parser.add_argument('-s', '--samples',
                        help='sample list with barcodes. \
                        If not provided then entire folder provided for --input is analyzed',
                        required=False)
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

    # demultiplex
    subparsers = parser.add_subparsers(dest="command", help='sub-command help')
    # Demultiplex
    sub_parser_demultiplex = subparsers.add_parser('demultiplex', help="Demultiplex reads into barcode samples")
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
    )

    # Genomeprep
    sub_parser_genome_prep = subparsers.add_parser('genomeprep', help="Prepare genome files from nucleotide sequences")
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
    )
    sub_parser_genome_prep.add_argument(
        "-gff", "--gff", help="generate GFF3 file", action="store_true", required=False
    )

    # Snakemake
    sub_parser_snakemake = subparsers.add_parser('snakemake', help="Run pyinseq using snakemake workflow")

    sub_parser_snakemake.add_argument('--get-default-config', action='store_true', help="")
    sub_parser_snakemake.add_argument('--config-format', default='yaml',
                       help="Format for the configuration file. Options are: 'json' or 'yaml'")
    sub_parser_snakemake.add_argument('-c', '--config', default=False)
    sub_parser_snakemake.add_argument('-T', '--threads', default=os.cpu_count(), type=int, help="Number of threads for running pyinseq")
    sub_parser_snakemake.add_argument('--additional-params',
                       help="Additional params passed to SNAKEMAKE. Make sure they are correct cause "
                            "I will not be checking....",
                       nargs='...',
                       default=[],
                       type=str,
                       )
    sub_parser_snakemake.add_argument(
        "-i", "--input", help="input Illumina reads file or folder", required=False
    )
    sub_parser_snakemake.add_argument(
        "-s", "--samples", help="sample list with barcodes", required=False
    )
    sub_parser_snakemake.add_argument(
        "-e",
        "--experiment",
        help="experiment name (no spaces or special characters)",
        required=False
    )
    sub_parser_snakemake.add_argument(
        "-g",
        "--genome",
        help="genome in GenBank format (one concatenated file for multiple contigs/chromosomes)",
        required=False,
    )
    sub_parser_snakemake.add_argument(
        "-d", "--disruption", help="fraction of gene disrupted (0.0 - 1.0)", default=1.0
    )
    sub_parser_snakemake.add_argument(
        "--min_count", help="Minimum number of reads per insertion site", default=3
    )
    sub_parser_snakemake.add_argument(
        "--max_ratio",
        help="Maximum ratio of left:right or right:left reads per insertion site",
        default=10,
    )
    return parser.parse_args()

