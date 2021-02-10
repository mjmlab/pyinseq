#!/usr/bin/env/python3
"""

Class and functions for running pyinseq using snakemake

"""
import subprocess as sub
from pathlib import Path
# Module imports
import pyinseq.workflows as w
from pyinseq.logger import logger
from pyinseq.runner import Settings
from pyinseq.utils import get_config_dict, create_experiment_directories


class PyinseqWorkflow():
    """
    Holds attributes and functions for running Snakemake
    """
    def __init__(self, config_file, add_params=['']):
        self.config_file = config_file
        self.config_dict = get_config_dict(self.config_file)
        self.experiment = self.config_dict['experiment']
        self.threads = self.config_dict['threads']
        self.snakefile = w.get_workflow_snakefile_path('PyinseqSnakemake')
        self.output_dir = Path(f"results/{self.experiment}")
        self.additional_params = add_params
        self.cmd = [
            'snakemake',
            '--configfile',
            self.config_file,
           '-s',
           self.snakefile,
           '--cores',
            str(self.threads),
            '--use-conda'
        ]
        # Add additional parameters
        self.cmd.extend(add_params)
        # Use pyinseq utils to create experiment directories, solves logging issue
        self.settings = Settings(self.experiment)
        create_experiment_directories(self.settings, snakemake=True)

    def run(self) -> None:
        """ Execute snakemake workflow using subprocess """
        logger.info("Calling snakemake for pyinseq workflow\n ")
        logger.info("All printing from pyinseq will be supressed but stored in "
                    f"{self.settings.summary_log}")
        # Call bash command, suppress pyinseq output
        sub.run(self.cmd, stdout=sub.DEVNULL)
        return





