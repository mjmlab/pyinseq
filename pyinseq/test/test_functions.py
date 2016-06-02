#!/usr/bin/env python
import pytest
from pyinseq import utils

# pyinseq.utils

def test_filename_replace_spaces():
    assert utils.convert_to_filename('file name') == 'file_name'

def test_filename_leading_whitespace():
    assert utils.convert_to_filename('  filename') == 'filename'

def test_filename_trailing_whitespace():
    assert utils.convert_to_filename('filename  ') == 'filename'

'''
# pyinseq.pyinseq

def test_set_disruption():
    assert pyinseq.set_disruption(1.0) == 1.0
'''
