#!/usr/bin/python
"""

Ad-hoc functions for snakemake workflows

"""

from pathlib import Path


def get_workflow_snakefile_path(workflow):
    """ Returns absolute path for Snakefile """
    return Path(__file__).parent.joinpath(f"{workflow}/Snakefile")
