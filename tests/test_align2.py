import pytest
from pydfinder import align2

def alignment_result():
	retun align2.AlignmentResult(...)

def test_is_valid_alignment(alignment_result):
	assert align2.is_valid_alignment(alignment_result, 0.05, 12)

def test_fitting_align():
	pass