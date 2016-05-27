#!/usr/bin/env python
import pytest
from pyinseq import utils

# pyinseq.utils.convert_to_filename

def test_filename_replace_spaces():
    assert utils.convert_to_filename('file name') == 'file_name'

def test_filename_leading_whitespace():
    assert utils.convert_to_filename('  filename') == 'filename'

def test_filename_trailing_whitespace():
    assert utils.convert_to_filename('filename  ') == 'filename'
