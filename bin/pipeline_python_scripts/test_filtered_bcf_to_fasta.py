from unittest import TestCase
from filtered_bcf_to_fasta import *
import tempfile

class BcfToFasta(TestCase):
    def test_calculate_reference_length(self):
        self.assertEqual(12, calculate_reference_length('test_reference.fas'))

    def test_calculate_gaps_to_add(self):
        self.assertEqual('-----', calculate_gaps_to_add(15,20))

    def test_filtered_bcf_to_fasta(self):
        self.assertEqual('---NacTN----', filtered_bcf_to_fasta('test_filtered.bcf', 12))

    def test_write_sequence(self):
        temp_filepath = tempfile.mkstemp()[1]
        sequence = 'GATCGATCGATC'
        write_sequence(temp_filepath, 'seq1', sequence)
        with open(temp_filepath) as input:
            lines = input.readlines()
            self.assertEqual('>seq1\n', lines[0])
            self.assertEqual('GATCGATCGATC\n', lines[1])