#!/usr/bin/env python
import unittest
from parse_csv import read_data

class OutputTests(unittest.TestCase):
    """These tests only work once the example01 data has been run."""

    def setUp(self):
        self.data = 'results/example01/summary_gene_table.txt'

    def test_csv_read_headers(self):
        self.assertEqual(
            read_data(self.data)[0],
            ['Contig', 'Start', 'End', 'Strand', 'Length', 'PID', 'Gene',
                'Synonym', 'Code', 'COG', 'Product', 'E001_01', 'E001_02'])

    def test_csv_row_with_transposon_hits(self):
        self.assertEqual(read_data(self.data)[12][6], 'gyrB')
        self.assertEqual(float(read_data(self.data)[12][11]), 1E5)
        self.assertEqual(float(read_data(self.data)[12][12]), 1E5)

    def test_csv_row_without_transposon_hits(self):
        self.assertEqual(read_data(self.data)[19][6], 'tusA')
        self.assertEqual(float(read_data(self.data)[19][11]), 0)
        self.assertEqual(float(read_data(self.data)[19][12]), 0)

    def test_csv_different_barcodes(self):
        self.assertEqual(read_data(self.data)[34][6], 'thiE')
        self.assertEqual(float(read_data(self.data)[34][11]), 5E4)
        self.assertEqual(float(read_data(self.data)[34][12]), 1E5)

if __name__ == '__main__':
    unittest.main()
