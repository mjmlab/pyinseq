#!/usr/bin/env python3
import csv
from collections import OrderedDict
from pyinseq.runner import Settings, tab_delimited_samples_to_dict
from pyinseq import utils
import numpy as np
import pandas as pd

# pyinseq.utils

def test_filename_replace_spaces():
    assert utils.convert_to_filename('file name') == 'file_name'

def test_filename_leading_whitespace():
    assert utils.convert_to_filename('  filename') == 'filename'

def test_filename_trailing_whitespace():
    assert utils.convert_to_filename('filename  ') == 'filename'

def test_filename_funny_characters():
    assert utils.convert_to_filename('filename*&*&*$$') == 'filename'
    assert utils.convert_to_filename('file*&*&*$$name') == 'filename'


# pyinseq.pyinseq

def test_class_Settings_init():
    x = Settings('example')
    assert x.path == 'results/example/'
    assert x.genome_path == 'results/example/genome_lookup/'
    assert x.raw_path == 'results/example/raw_data/'
    assert x.samples_yaml == 'results/example/samples.yml'
    assert x.summary_yaml == 'results/example/summary.yml'
    assert x.keepall == False
    assert x.barcode_length == 4

def test_tab_delimited_samples_to_dict_trailing_newline():
    s = 'pyinseq/tests/data/additional/sample01_01.txt'
    assert tab_delimited_samples_to_dict(s) == \
        OrderedDict([('sample_1', {'barcode': 'AAAA'}), ('sample_2', {'barcode': 'TTTT'})])

def test_tab_delimited_samples_to_dict_no_trailing_newline():
    s = 'pyinseq/tests/data/additional/sample01_02.txt'
    assert tab_delimited_samples_to_dict(s) == \
        OrderedDict([('sample_1', {'barcode': 'AAAA'}), ('sample_2', {'barcode': 'TTTT'})])
