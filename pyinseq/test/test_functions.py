#!/usr/bin/env python
import unittest
from pyinseq import utils

class FunctionTests(unittest.TestCase):

    # pyinseq.utils.convert_to_filename

    def test_filename_replace_spaces(self):
        self.assertEqual(utils.convert_to_filename('file name'), 'file_name')

    def test_filename_leading_whitespace(self):
        self.assertEqual(utils.convert_to_filename('  filename'), 'filename')

    def test_filename_trailing_whitespace(self):
        self.assertEqual(utils.convert_to_filename('filename  '), 'filename')


if __name__ == '__main__':
    unittest.main()
