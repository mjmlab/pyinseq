#!/usr/bin/env/python3

"""

Settings for running baseline pyinseq

"""

import os
import yaml
from pathlib import Path
from pyinseq.logger import pyinseq_logger as logger
from pyinseq.workflows import get_workflow_snakefile_path
from pyinseq.utils import \
    convert_to_filename, \
    tab_delimited_samples_to_dict, \
    get_config_dict


class Settings:
    """Instantiate to set up settings for the experiment."""
    def __init__(self, command, config_file, add_params=[]):
        # command options are: ['pyinseq', 'demultiplex']
        self.command = command if command else 'pyinseq'
        # Attributes for all workflows
        self.config_file = config_file
        config_dict = get_config_dict(self.config_file)
        self.experiment = convert_to_filename(config_dict['experiment'])
        self.threads = config_dict['threads']
        self.output_dir = Path(f"results/{self.experiment}")
        self.path = f"results/{self.experiment}/"
        self.log = f"{self.path}log.txt"
        self.summary_log = f"{self.path}summary_log.txt"

        if self.command == 'pyinseq':
            self.reads = config_dict['input']
            self.samples = config_dict['samples']
            self.samples_txt = f"{self.path}samples.txt"
            self.samples_dict = tab_delimited_samples_to_dict(self.samples)
            self.samples_info_yml = f"{self.path}samples_info.yml"
            self.raw_path = f"{self.path}raw_data/"
            self.reference_genome = config_dict['genome']
            self.genome_path = f"{self.path}genome_lookup/"
            self.map_to_genome = True
            self.gff3 = config_dict['gff3']
            self.summary_table = f"{self.path}summary_gene_table.txt"
            # Pyinseq optional args
            self.disruption = set_disruption(config_dict['disruption'])
            self.barcode_length = set_barcode_length(config_dict['barcode_length'])
            self.transposon_seq = set_transposon_seq(config_dict['transposon_seq'])
            # counts at one transposon site for it to qualify
            self.min_counts = set_min_count(config_dict['min_count'])
            # max ratio of left/right sites for it to qualify
            self.max_ratio = set_max_ratio(config_dict['max_ratio'])

        elif self.command == 'demultiplex':
            self.reads = config_dict['input']
            self.samples_txt = f"{self.path}samples.txt"
            self.samples = config_dict['samples']
            self.samples_dict = tab_delimited_samples_to_dict(self.samples)
            self.samples_info_yml = f"{self.path}samples_info.yml"
            self.raw_path = f"{self.path}raw_data/"
            self.barcode_length = set_barcode_length(config_dict['barcode_length'])
            self.transposon_seq = set_transposon_seq(config_dict['transposon_seq'])

        elif self.command == 'genomeprep':
            self.reference_genome = config_dict['genome']
            self.genome_path = f"{self.path}genome_lookup/"
            self.gff3 = config_dict['gff3']

        self.settings_pickle = self.output_dir.joinpath('settings.pickle')
        # organism reference files called 'genome.fna' etc
        self.organism = "genome"

        # Set shell command
        self.snakefile = get_workflow_snakefile_path(self.command)
        self.snakemake_cmd = [
            'snakemake',
            '--config',
            f'experiment={self.experiment}',
            '-s',
            str(self.snakefile),
            '--cores',
            str(self.threads - 1),
        ]
        # Add additional parameters, if any...
        self.snakemake_cmd.extend(add_params)
        self.snakemake_cmd = ' '.join(self.snakemake_cmd)
        # Create experiment directories, solves logging issue
        self.process_reads = True if command in ['pyinseq', 'demultiplex'] else False
        self.process_sample_list = True if command in ['pyinseq', 'demultiplex'] else False
        self.parse_genebank = True if command in ['pyinseq', 'genomeprep'] else False
        self.generate_bowtie_index = True if command in ['pyinseq', 'genomeprep'] else False
        self.write_trimmed_reads = not config_dict['notrim'] if command == 'demultiplex' else True

    def __repr__(self):
        # Print each variable on a separate line
        return "\n".join("%s: %s" % item for item in vars(self).items())

    def clean_up(self):
        """ Cleans up settings pickle file"""
        try:
            os.remove(self.settings_pickle)
        except OSError:
            pass

    def dump_sample_dict_to_yml(self, updated_dict=None):
        """ Updates sample_info.yml file with attributes for sample """
        if not updated_dict:
            # Dump yml file for info
            with open(self.samples_info_yml, 'w') as f:
                yaml.dump(self.samples_dict, f)
            return
        # Load current state of sample_dict
        with open(self.samples_info_yml, 'r+') as f:
            current_dict = yaml.load(f, Loader=yaml.FullLoader)
        # Updated the old sample with new information
        with open(self.samples_info_yml, "w") as fo:
            dump_dict = current_dict.copy()
            for sample, sub_d in updated_dict.items():
                for k in sub_d:
                    if k not in current_dict[sample].keys():
                        dump_dict[sample][k] = updated_dict[sample][k]
            yaml.dump(dump_dict, fo)
        return


def set_min_count(min_count=3):
    """Check min_count is positive, otherwise set default."""
    if min_count < 0:
        logger.error(
            f"Min_count value provided ({min_count}) is not positive; proceeding with default value of 3"
        )
        return 3
    return min_count


def set_max_ratio(max_ratio=10):
    """Check that max_ratio is positive, otherwise set default."""
    if max_ratio < 0:
        logger.error(
            f"Max_ratio value provided ({max_ratio}) is not positive; proceeding with default value of 10"
        )
        return 10
    return max_ratio


def set_disruption(d=0.9):
    """Check that gene disruption is 0.0 to 1.0; otherwise set to 0.9."""
    if d < 0.0 or d > 1.0:
        logger.error(
            f"Disruption value provided ({d}) is not in range 0.0 to 0.9; proceeding with default value of 0.9"
        )
        return 0.9
    return d


def set_barcode_length(barcode_length=4):
    """Check that min_count and max_ratio are positive, otherwise set default."""
    if barcode_length > 16:
        logger.error(
            f"Barcode length {barcode_length} is larger than 16bp which is unsupported by pyinseq. "
            f"Setting length to 4bp default"
        )
        return 4
    if barcode_length < 4:
        logger.error(
            f"Barcode length {barcode_length} is smaller than 4bp which is unsupported by pyinseq. "
            f"Setting length to 4bp default"
        )
        return 4
    return barcode_length


def set_transposon_seq(transposon_seq='ACAGGTTG'):
    """Check that min_count and max_ratio are positive, otherwise set default."""
    if not transposon_seq.isalpha():
        logger.error(
            f"Transposon sequence provided contains non-letter characters. "
            f"Setting default to ACAGGTTG"
        )
        return 'ACAGGTTG'
    return transposon_seq

