#!/usr/bin/env python3
from .test_utils import runscript, datadir
import filecmp
import pytest
import multiprocessing as mp # TEST
import os # TEST
import pickle # TEST

def test_pyinseq_script_no_args(datadir,tmpdir):

    args = []
    status, out, err = runscript('pyinseq', args, directory=str(tmpdir))
    assert out[0:5] == 'usage'
    print(status,err)

def test_pyinseq_script(datadir, tmpdir):
    '''
    input_fn = datadir('input/example01.fastq')
    sample_fn = datadir('input/example01.txt')
    gb_fn = datadir('input/ES114v2.gb')
    output_name = 'example_pyinseq'
    expected_output = datadir('output_pyinseq')
    output_dir = tmpdir.join('results/example_pyinseq')

    args = ['-i', input_fn, '-s', sample_fn, '-g', gb_fn, '-e', output_name]
    status, out, err = runscript('pyinseq', args, directory=str(tmpdir))
    print(status, out, err)
    assert status  

    with open("pickle_objects","wb") as pick:
        pickle.dump([expected_output,output_dir],pick)
    '''
    
    with open("pickle_objects",'rb') as pick: # Records local variables as dictionary for efficient debugging
        pickle_dict = pickle.load(pick)

    expected_output = pickle_dict[0]; output_dir = pickle_dict[1]

    dcmp = filecmp.dircmp(expected_output,
                          output_dir,
                          ignore=['E001_01_bowtie.txt', 'E001_02_bowtie.txt','.DS_Store','log.txt'])

    assert 'log.txt' in os.listdir(output_dir) # check that log file is created from pyinseq

    # checks that files are same in both directories
    assert not dcmp.left_only and not dcmp.right_only
    
    # check files to see if content differs
    assert not dcmp.diff_files 
    assert not dcmp.funny_files # Check for files that cannot be compared

    # because subdirs is a dict keyed by subdir name
    # with dircmp objects as values
    for subdcmp in dcmp.subdirs.values():
        assert not subdcmp.diff_files
        assert not subdcmp.funny_files
        assert not subdcmp.left_only and not subdcmp.right_only
    

def test_pyinseq_demultiplex_script(datadir, tmpdir):

    input_fn = datadir('input/example01.fastq')
    sample_fn = datadir('input/example01.txt')
    output_name = 'example_demultiplex'
    expected_output = datadir('output_demultiplex')
    output_dir = tmpdir.join('results/example_demultiplex')

    args = ['demultiplex', '-i', input_fn, '-s', sample_fn, '-e', output_name]
    status, out, err = runscript('pyinseq', args, directory=str(tmpdir))

    assert status == 0

    dcmp = filecmp.dircmp(expected_output,
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

    dcmp = filecmp.dircmp(expected_output,
                          str(output_dir))
    assert dcmp.diff_files == []

    # because subdirs is a dict keyed by dir name
    # with dircmp object values
    for subdcmp in dcmp.subdirs.values():
        assert subdcmp.diff_files == []


def test_pyinseq_genomeprep_script(datadir, tmpdir):

    output_name = 'example_genomeprep'
    expected_output = datadir('output_genomeprep')
    output_dir = tmpdir.join('results/example_genomeprep')
    gb_fn = datadir('input/ES114v2.gb')

    args = ['genomeprep', '-e', output_name, '-g', gb_fn]
    status, out, err = runscript('pyinseq', args, directory=str(tmpdir))

    assert status == 0

    dcmp = filecmp.dircmp(expected_output,
                          str(output_dir))
    assert dcmp.diff_files == []
    # because subdirs is a dict keyed by dir name
    # with dircmp object values
    for subdcmp in dcmp.subdirs.values():
        assert subdcmp.diff_files == []