#!/usr/bin/env python3
from .test_utils import runscript, datadir
import filecmp
import pytest

def test_pyinseq_script(datadir, tmpdir):

    input_fn = datadir('input/example01.fastq')
    sample_fn = datadir('input/example01.txt')
    gb_fn = datadir('input/ES114v2.gb')
    output_name = 'example'
    expected_output = datadir('output')
    output_dir = tmpdir.join('results/example')

    args = ['-i', input_fn, '-s', sample_fn, '-g', gb_fn, '-e', output_name]
    status, out, err = runscript('pyinseq', args, directory=str(tmpdir))

    assert status == 0

    dcmp = filecmp.dircmp(datadir('output'),
                          str(output_dir),
                          ignore=['E001_01_bowtie.txt', 'E001_02_bowtie.txt'])
    assert dcmp.diff_files == []
    # because subdirs is a dict keyed by dir name
    # with dircmp object values
    for subdcmp in dcmp.subdirs.values():
        assert subdcmp.diff_files == []
