import pytest
from pydfinder import align2
from pydfinder import align

# def alignment_result():
#   retun align2.AlignmentResult(...)

# def test_is_valid_alignment(alignment_result):
#   assert align2.is_valid_alignment(alignment_result, 0.05, 12)

def test_fitting_align():
    v = "GTAGGCTTAAGGTTA"
    w = "TAGATA"
    res = align2.AlignmentResult(1, 9, 0, 6,"TAGGCTTA", "TAGA-T-A", 2)
    assert align2.fitting_alignment(v, w, align.ScoreParam(1, -1, 1)) == res

def test_good_call_from_5e_soft_clipping():
    read = align2.ForwardFiveEndSoftClipping()
    call = Call()
    assert read.find_call(target_seq, aligner) == call

def test_null_call_from_5e_soft_clipping():
    pass
