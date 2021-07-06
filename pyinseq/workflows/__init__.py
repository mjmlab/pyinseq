#!/usr/bin/python

"""

Ad-hoc functions for snakemake workflows

"""

from pathlib import Path


def get_workflow_snakefile_path(command):
    """Returns absolute path for Snakefile"""
    workflows = {
        "pyinseq": Path(__file__).parent.joinpath(f"PyinseqWorkflow/Snakefile"),
        "demultiplex": Path(__file__).parent.joinpath(f"DemultiplexWorkflow/Snakefile"),
        "genomeprep": Path(__file__).parent.joinpath(f"GenomeprepWorkflow/Snakefile"),
    }
    return workflows[command]
