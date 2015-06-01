import difflib

def overlaps(s1, s2, min_len, min_sim):
	l = len(s1)
	assert(min_len <= l)
	for i in range(min_len, l+1):
		sm = difflib.SequenceMatcher(None, s1[-i:], s2[:i])
		if sm.quick_ratio > min_sim:
			return True
		sm2 = difflib.SequenceMatcher(None, s1[:i], s2[-i:])
		if sm2.quick_ratio > min_sim:
			return True
	return False
