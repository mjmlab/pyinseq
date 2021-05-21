#!/usr/bin/env python3

"""

Tests for individual functions in pyinseq module

"""

import filecmp
import os
from collections import OrderedDict

# Setup logger before importing modules
from pyinseq.logger import pyinseq_logger

pyinseq_logger.setup_logger()
# Module imports
from pyinseq import utils
from pyinseq.settings import Settings
from .test_utils import cd, runscript, datadir, load_settings


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
    setting = load_settings

    assert setting.path == "results/test_pyinseq/"
    assert setting.genome_path == "results/test_pyinseq/genome_lookup/"
    assert setting.raw_path == "results/test_pyinseq/raw_data/"
    assert setting.samples_txt == "results/test_pyinseq/samples.txt"
    assert setting.log == "results/test_pyinseq/log.txt"
    assert setting.summary_log == "results/test_pyinseq/summary_log.txt"
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
    settings = load_settings
    expected_output = datadir("output_pyinseq")
    # In the context of pytest directory
    with cd(str(tmpdir)):
        # Manually modify settings object
        from pathlib import Path

        settings.output_dir = Path(str(tmpdir) + "results/test_pyinseq")
        utils.create_experiment_directories(settings)
        # Check /results/experiment exist
        assert os.path.exists(settings.output_dir)
        assert filecmp.dircmp(expected_output, settings.output_dir).common_dirs == [
            "genome_lookup",
            "raw_data",
        ]
