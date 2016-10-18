#!/usr/bin/env python
import pytest
from .parse_csv import read_data

# Set up

example01summary = read_data('results/example01/summary_gene_table.txt')

def test_csv_read_headers():
    # data = read_data('results/example01/summary_gene_table.txt')
    assert example01summary[0] == ['Contig', 'Start', 'End', 'Strand', 'Length',
                                   'PID', 'Gene', 'Synonym', 'Code', 'COG',
                                   'Product', 'E001_01', 'E001_02']

def test_csv_row_with_transposon_hits():
    assert example01summary[12][6] == 'gyrB'
    assert float(example01summary[12][11]) == 1E5
    assert float(example01summary[12][12]) == 1E5

def test_csv_row_without_transposon_hits():
    assert example01summary[19][6] == 'tusA'
    assert float(example01summary[19][11]) == 0
    assert float(example01summary[19][12]) == 0

def test_csv_with_different_data_in_different_barcoded_samples():
    assert example01summary[34][6] == 'thiE'
    assert float(example01summary[34][11]) == 5E4
    assert float(example01summary[34][12]) == 1E5
