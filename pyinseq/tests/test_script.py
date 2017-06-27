#!/usr/bin/env python3
from .test_utils import runscript, datadir
import filecmp
import pytest


def test_pyinseq_script(datadir, tmpdir):

    input_fn = datadir('input/example01.fastq')
    sample_fn = datadir('input/example01.txt')
    gb_fn = datadir('input/ES114v2.gb')
    output_name = 'example_pyinseq'
    expected_output = datadir('output_pyinseq')
    output_dir = tmpdir.join('results/example_pyinseq')

    args = ['-i', input_fn, '-s', sample_fn, '-g', gb_fn, '-e', output_name]
    status, out, err = runscript('pyinseq', args, directory=str(tmpdir))

    assert status == 0

    dcmp = filecmp.dircmp(datadir('output_pyinseq'),
                          str(output_dir),
                          ignore=['E001_01_bowtie.txt', 'E001_02_bowtie.txt'])
    assert dcmp.diff_files == []
    # because subdirs is a dict keyed by dir name
    # with dircmp object values
    for subdcmp in dcmp.subdirs.values():
        assert subdcmp.diff_files == []


def test_pyinseq_demultiplex_script(datadir, tmpdir):

    input_fn = datadir('input/example01.fastq')
    sample_fn = datadir('input/example01.txt')
    output_name = 'example_demultiplex'
    expected_output = datadir('output_demultiplex')
    output_dir = tmpdir.join('results/example_demultiplex')

    args = ['demultiplex', '-i', input_fn, '-s', sample_fn, '-e', output_name]
    status, out, err = runscript('pyinseq', args, directory=str(tmpdir))

    assert status == 0

    dcmp = filecmp.dircmp(datadir('output_demultiplex'),
                          str(output_dir))
    assert dcmp.diff_files == []
    # because subdirs is a dict keyed by dir name
    # with dircmp object values
    for subdcmp in dcmp.subdirs.values():
        assert subdcmp.diff_files == []


def test_pyinseq_demultiplex_notrim_script(datadir, tmpdir):

    input_fn = datadir('input/example01.fastq')
    sample_fn = datadir('input/example01.txt')
    output_name = 'example_demultiplex_notrim'
    expected_output = datadir('output_demultiplex_notrim')
    output_dir = tmpdir.join('results/example_demultiplex_notrim')

    args = ['demultiplex', '-i', input_fn, '-s', sample_fn, '-e', output_name, '--notrim']
    status, out, err = runscript('pyinseq', args, directory=str(tmpdir))

    assert status == 0

    dcmp = filecmp.dircmp(datadir('output_demultiplex_notrim'),
                          str(output_dir))
    assert dcmp.diff_files == []
    # because subdirs is a dict keyed by dir name
    # with dircmp object values
    for subdcmp in dcmp.subdirs.values():
        assert subdcmp.diff_files == []
