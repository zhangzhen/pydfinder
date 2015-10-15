import unittest
from pydfinder import align

class AlignerTest(unittest.TestCase):

	maxDiff = None

	def setUp(self):
		v = "ACGTTACT"
		w = "ACGGGGACT"
		self.aligner = align.Aligner(v, w)

	def test_calc_score_matrix(self):
		S = [[0,0,0,0,0,0,0,0,0,0],
		[0,1,0,0,0,0,0,1,0,0],
		[0,0,2,0,0,0,0,0,2,0],
		[0,0,0,3,1,1,1,0,0,1],
		[0,0,0,0,2,0,0,0,0,1],
		[0,0,0,0,0,1,0,0,0,1],
		[0,1,0,0,0,0,0,1,0,0],
		[0,0,2,0,0,0,0,0,2,0],
		[0,0,0,1,0,0,0,0,0,3]]

		self.aligner._calc_score_matrix(align.ScoreParam(1, -1, 2, 4))
		self.assertEqual(S, self.aligner.S)

	def test_find_max_pos(self):
		self.aligner._calc_score_matrix(align.ScoreParam(1, -1, 2, 4))
		self.aligner._calc_max_matrix()
		self.assertEqual((3,3), self.aligner._find_max_pos())
		