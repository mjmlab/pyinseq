#!/usr/bin/env python3

"""

Tests for individual functions in pyinseq module

"""

import os
import sys
import filecmp
from pathlib import Path
from collections import OrderedDict

# Setup logger before importing modules
from pyinseq.logger import pyinseq_logger

pyinseq_logger.setup_logger()
# Module imports
from pyinseq import utils
from pyinseq.settings import Settings
from pyinseq.parsers import get_args

# Test imports
from .test_utils import cd, runscript, datadir, load_settings, get_dump


def test_filename_replace_spaces():
    assert utils.convert_to_filename("file name") == "file_name"


def test_filename_leading_whitespace():
    assert utils.convert_to_filename("  filename") == "filename"


def test_filename_trailing_whitespace():
    assert utils.convert_to_filename("filename  ") == "filename"


def test_filename_funny_characters():
    assert utils.convert_to_filename("filename*&*&*$$") == "filename"
    assert utils.convert_to_filename("file*&*&*$$name") == "filename"


def test_pyinseq_Settings(datadir, load_settings):
    setting = load_settings("pyinseq")

    assert setting.path == Path("results/test_pyinseq/")
    assert setting.genome_path == Path("results/test_pyinseq/genome_lookup/")
    assert setting.raw_path == Path("results/test_pyinseq/raw_data/")
    assert setting.samples_txt == Path("results/test_pyinseq/samples.txt")
    assert setting.log == Path("results/test_pyinseq/log.txt")
    assert setting.summary_log == Path("results/test_pyinseq/summary_log.txt")
    assert setting.barcode_length == 4


def test_tab_delimited_samples_to_dict_trailing_newline():
    s = "pyinseq/tests/data/additional/sample01_01.txt"
    assert utils.tab_delimited_samples_to_dict(s) == OrderedDict(
        [("sample_1", {"barcode": "AAAA"}), ("sample_2", {"barcode": "TTTT"})]
    )


def test_tab_delimited_samples_to_dict_no_trailing_newline():
    s = "pyinseq/tests/data/additional/sample01_02.txt"
    assert utils.tab_delimited_samples_to_dict(s) == OrderedDict(
        [("sample_1", {"barcode": "AAAA"}), ("sample_2", {"barcode": "TTTT"})]
    )


def test_create_experiment_directories(datadir, tmpdir, load_settings):
    """Test for function that creates the directory structure of pyinseq"""
    settings = load_settings("pyinseq")
    expected_output = datadir("output_pyinseq")
    # In the context of pytest directory
    with cd(str(tmpdir)):
        # Manually modify settings object
        from pathlib import Path

        settings.path = Path(tmpdir + "results/test_pyinseq")
        utils.create_experiment_directories(settings)
        # Check /results/experiment exist
        assert os.path.exists(settings.path)
        assert filecmp.dircmp(expected_output, settings.path).common_dirs == [
            "genome_lookup",
            "raw_data",
        ]


def test_write_config_file(datadir, tmpdir):
    """Test for writing configuration file"""
    input_fn = "../data/input/example01.fastq"
    sample_fn = "../data/input/example01.txt"
    gb_fn = "../data/input/ES114v2.gb"
    output_name = "test_pyinseq"
    expected_config = datadir("input/test-pyinseq-config.yml")
    output_config = tmpdir.join("test_pyinseq-config.yml")

    # Modify sys.argv
    sys.argv = [
        "pyinseq",
        "-i",
        input_fn,
        "-s",
        sample_fn,
        "-g",
        gb_fn,
        "-e",
        output_name,
        "-d",
        "0.9",
        "--threads",
        "2",
        "--snakemake_params",
        "--use-conda",
    ]
    # Get the Namespace object from argparse
    parser, args = get_args()
    with cd(str(tmpdir)):
        utils.write_config_file(args)
        # Check files
        assert filecmp.cmp(expected_config, output_config)
