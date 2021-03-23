#!/usr/bin/python
"""

Ad-hoc functions for snakemake workflows

"""

import subprocess as sub
from pathlib import Path
# Module imports
from pyinseq.settings import Settings
from pyinseq.utils import get_config_dict, create_experiment_directories


class PyinseqWorkflow:
    """

    Class with base attributes for pyinseq workflows

    """
    def __init__(self, command, config_file, add_params=['']):
        self.command = command if command else 'pyinseq'
        self.snakefile = get_workflow_snakefile_path(self.command)
        self.config_file = config_file
        self.config_dict = get_config_dict(self.config_file)
        self.experiment = self.config_dict['experiment']
        self.threads = self.config_dict['threads']
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
        ]
        # Add additional parameters, if any...
        self.cmd.extend(add_params)
        # Create experiment directories, solves logging issue
        self.settings = Settings(self.experiment)
        create_experiment_directories(self.settings)

    def execute(self) -> None:
        """ Execute snakemake workflow using subprocess """
        # Call bash command, suppress pyinseq output
        sub.run(self.cmd)
        return


def get_workflow_snakefile_path(command):
    """ Returns absolute path for Snakefile """
    workflows = {
        'pyinseq': Path(__file__).parent.joinpath(f"PyinseqWorkflow/Snakefile"),
        'demultiplex': Path(__file__).parent.joinpath(f"DemultiplexWorkflow/Snakefile"),
        'genomeprep': Path(__file__).parent.joinpath(f"GenomeprepWorkflow/Snakefile"),
    }
    return workflows[command]
