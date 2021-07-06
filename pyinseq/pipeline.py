#!/usr/bin/env python3

"""

Main script for running the pyinseq package

"""

import yaml

# Module imports
from pyinseq.utils import copy_file
from pyinseq.analyze import t_fifty
from pyinseq.settings import Settings
from pyinseq.logger import pyinseq_logger
from pyinseq.map_reads import bowtie_map, parse_bowtie
from pyinseq.process_mapping import map_sites, map_genes, build_gene_table

logger = pyinseq_logger.logger


def pipeline_mapping(settings, samples_dict):
    """Aggregate bowtie output, map to genes in the feature table, and aggregate samples."""
    # Dictionary of each sample's cpm by gene
    gene_mappings = {}
    mapping_data = {}
    for sample in samples_dict:
        with cd(settings.genome_path):
            # Paths are relative to the genome_lookup directory
            # from where bowtie is called
            bowtie_in = "../" + sample + "_trimmed.fastq"
            bowtie_out = "../" + sample + "_bowtie.txt"
            # map to bowtie and produce the output file
            logger.info(f"Sample {sample}: map reads with bowtie")
            bowtie_msg_out = bowtie_map(settings.organism, bowtie_in, bowtie_out)
            logger.info(bowtie_msg_out)
            # store bowtie data for each sample in dictionary
            mapping_data[sample] = {"bowtie_results": [], "insertion_sites": []}
            mapping_data[sample]["bowtie_results"] = parse_bowtie(bowtie_msg_out)
        # Map each bowtie result to the chromosome
        insertions = len(map_sites(sample, settings))
        mapping_data[sample]["insertion_sites"] = insertions
        # Add gene-level results for the sample to geneMappings
        # Filtered on gene fraction disrupted as specified by -d flag
        gene_mappings[sample] = map_genes(sample, settings)
        # if not settings.keep_all:
        #    # Delete trimmed fastq file, bowtie mapping file after writing mapping results
        #    os.remove(s["trimmedPath"])
        #    os.remove("results/{Settings.experiment}/{bowtieOutputFile}"
        t_fifty_result = t_fifty(sample, settings)
        logger.info(f"T50 result for {sample}: {t_fifty_result}")
    build_gene_table(
        settings.organism, samples_dict, gene_mappings, settings.experiment
    )


def pipeline_summarize(samples_dict: dict, settings: Settings):
    """Summary of INSeq run."""
    # Write summary log with more data
    copy_file(settings.samples, settings.samples_txt)
    logger.info(f"Print log path: {settings.log}")
    pyinseq_logger.logger_io.write(f"Print log path: {settings.log}\n")
    logger.info(f"Print settings:\n{settings}\n")
    pyinseq_logger.logger_io.write(f"Print settings:\n{settings}\n")
    with open(settings.samples_info_yml, "r") as f:
        d = yaml.load(f, Loader=yaml.FullLoader)
    logger.info(f"Print samples detail\n{yaml.dump(d, default_flow_style=False)}\n")


def pipeline_analysis(samples_dict: dict, settings: Settings) -> None:
    """Pipeline for analyzing pyinseq results"""
    # T50 calculation
    T50_dict = dict()
    for sample in samples_dict:
        T50_dict[sample] = t_fifty(sample, settings)
        logger.info(f"T50 {sample}: {T50_dict[sample]}")
        # plot_insertions(sample, settings)
    return T50_dict


if __name__ == "__main__":
    pass
