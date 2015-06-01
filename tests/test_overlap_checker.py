import unittest
from dfinder2 import overlap_checker

class TestOverlapChecker(unittest.TestCase):

	def test_plus_overlaps(self):
		s1 = "AGGGATGGAATAAAAGCCCTTTTTGCTTCCGTTACAGAGTCCTGCAGTTCTGTTTGAATATTAATTAATTCCCGACCTGATCCTCATCTTCTCAAGGAACA"
		s2 = "GCAGTTCTGTTTGAATATTAATCAATTCCCGACCTGATCCTCATCTTCTCAAGGAACAATTGCAGTGAATCCTCAGTTGCCTACGTGGTTAATGTGGAAGA"
		min_len = 37
		min_sim = 0.95
		self.assertTrue(overlap_checker.overlaps(s1, s2, min_len, min_sim))

	def test_minus_overlaps(self):
		s1 = "AGGGATGGAATAAAAGCCCTTTTTGCTTCCGTTACAGAGTCCTGCAGTTCTGTTTGAATATTAATTAATTCCCGACCTGATCCTCATCTTCTCAAGGAACA"
		s2 = "GCAGTTCTGTTTGAATATTAATCAATTCCCGACCTGATCCTCATCTTCTCAAGGAACAATTGCAGTGAATCCTCAGTTGCCTACGTGGTTAATGTGGAAGA"
		min_len = 37
		min_sim = 0.95
		self.assertTrue(overlap_checker.overlaps(s2, s1, min_len, min_sim))
