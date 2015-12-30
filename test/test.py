import unittest
from parse_csv import read_data

class OutputTests(unittest.TestCase):

    def setUp(self):
        self.data = 'example01/summary_gene_table.txt'

    def test_csv_read_headers(self):
        self.assertEqual(
            read_data(self.data)[0],
            ['Contig', 'Start', 'End', 'Strand', 'Length', 'PID', 'Gene',
            'Synonym', 'Code', 'COG', 'Product', 'E001_01', 'E001_02']
            )

    #def test_csv_read_data_points(self):
    #    self.assertEqual(read_data(self.data)[1][7], '87')

if __name__ == '__main__':
    unittest.main()
