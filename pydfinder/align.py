DIAGONAL = 1
VERTICAL = 2
HORIZONTAL = 3

def make_matrix(sizex, sizey):
    """Creates a sizex by sizey matrix filled with zeros."""
    return [[0]*sizey for i in xrange(sizex)]

#=============================================================
# Alignment Parameters
#=============================================================

class ScoreParam(object):
    """Stores the parameters for an alignment scoring function"""
    def __init__(self, match, mismatch, gap, gap_start=0):
        self.gap_start = gap_start
        self.gap = gap
        self.match = match
        self.mismatch = mismatch

    def matchchar(self, a,b):
        """Return the score for aligning character a with b"""
        assert len(a) == len(b) == 1
        if a==b:
            return self.match
        else:
            return self.mismatch

    def __str__(self):
        return "match = %d; mismatch = %d; gap_start = %d; gap_extend = %d" % (
                self.match, self.mismatch, self.gap_start, self.gap
        )


class Aligner(object):
    """docstring for Aligner"""
    def __init__(self, v, w):
        self.v = v
        self.w = w
        self.size_v = len(v)
        self.size_w = len(w)

    def local_align(self, score_param):
        _calc_score_matrix(score_param)
        _calc_max_matrix()
        max_i, max_j = _find_max_pos()
        v_aligned, w_aligned = backtrace(max_i, max_j)
        return (max_i, max_j), v_aligned, w_aligned
    
    def _calc_score_matrix(self, score_param):
        self.S = make_matrix(self.size_v + 1, self.size_w + 1)
        self.S_backtrace = make_matrix(self.size_v + 1, self.size_w + 1)
        S_lower = make_matrix(self.size_v + 1, self.size_w + 1)
        S_upper = make_matrix(self.size_v + 1, self.size_w + 1)

        # Fill in the Score and Backtrace matrices.
        for i in xrange(1, self.size_v+1):
            for j in xrange(1, self.size_w+1):
                S_lower[i][j] = max([S_lower[i-1][j] - score_param.gap, self.S[i-1][j] - score_param.gap_start])
                S_upper[i][j] = max([S_upper[i][j-1] - score_param.gap, self.S[i][j-1] - score_param.gap_start])
                middle_scores = [0, self.S[i-1][j-1] + score_param.matchchar(self.v[i-1], self.w[j-1]), S_lower[i][j], S_upper[i][j]]
                self.S[i][j] = max(middle_scores)
                self.S_backtrace[i][j] = middle_scores.index(self.S[i][j])

    def _calc_max_matrix(self):
        self.M_backtrace = make_matrix(self.size_v + 1, self.size_w + 1)

        for i in xrange(1, self.size_v + 1):
            self.M_backtrace[i][0] = VERTICAL

        for i in xrange(1, self.size_w + 1):
            self.M_backtrace[0][i] = HORIZONTAL

        for i in xrange(1, self.size_v + 1):
            for j in xrange(1, self.size_w + 1):
                scores = [self.S[i][j], self.S[i-1][j], self.S[i][j-1]]
                self.S[i][j] = max(scores)
                if scores.index(self.S[i][j]) != 0:
                    self.M_backtrace[i][j] = scores.index(self.S[i][j]) + 1

    def _find_max_pos(self):
        max_score = self.S[self.size_v][self.size_w]

        for i in xrange(0, self.size_v + 1):
            for j in xrange(0, self.size_w + 1):
                if self.M_backtrace[i][j] == 0 and self.S[i][j] == max_score:
                    return (i, j)

        return None

    def _backtrace(self, i, j):

        insert_indel = lambda word, i: word[:i] + '-' + word[i:]

        v_aligned, w_aligned = self.v[:j], self.w[:i]

        while self.S_backtrace[i][j] != 0 and i*j != 0:
            if backtrace[i][j] == VERTICAL:
                i -= 1
                w_aligned = insert_indel(w_aligned, j)
            elif backtrace[i][j] == HORIZONTAL:
                j -= 1
                v_aligned = insert_indel(v_aligned, i)
            elif backtrace[i][j] == DIAGONAL:
                i -= 1
                j -= 1

        v_aligned = v_aligned[i:]
        w_aligned = w_aligned[j:]

        return (v_aligned, w_aligned)
