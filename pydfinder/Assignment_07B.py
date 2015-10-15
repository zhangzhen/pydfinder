#!/usr/bin/env python
'''
A solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
The associated textbook is Bioinformatics Algorithms: An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.
The course is run on Coursera and the assignments and textbook are hosted on Stepic

Problem Title: Fitting Alignment Problem
Assignment #: 07
Problem ID: B
URL: https://stepic.org/Bioinformatics-Algorithms-2/The-Changing-Faces-of-Sequence-Alignment-248/step/5
'''

import itertools
import operator

def hamming(str1, str2):
    assert len(str1) == len(str2), "length should be equal:\n{}\n{}".format(str1, str2)
    #ne = str.__ne__  ## this is surprisingly slow
    ne = operator.ne
    return sum(itertools.imap(ne, str1, str2))

def fitting_alignment(v,w):
    '''Returns the fitting alignment of strings v and w, along with the associated score.'''
    # Initialize the matrices.
    S = [[0 for j in xrange(len(w)+1)] for i in xrange(len(v)+1)]
    backtrack = [[0 for j in xrange(len(w)+1)] for i in xrange(len(v)+1)]

    for j in xrange(1, len(w)+1):
        S[0][j] = -3*j

    # Fill in the Score and Backtrack matrices.
    for i in xrange(1, len(v)+1):
        for j in xrange(1, len(w)+1):
            scores = [S[i-1][j] - 3, S[i][j-1] - 3, S[i-1][j-1] + [-1, 1][v[i-1] == w[j-1]]]
            S[i][j] = max(scores)
            backtrack[i][j] = scores.index(S[i][j])

    # Get the position of the highest scoring cell corresponding to the end of the shorter word w.
    j = len(w)
    i = max(enumerate([S[row][j] for row in xrange(len(w), len(v))]),key=lambda x: x[1])[0] + len(w)
    max_score = S[i][j]
    w_e = j
    v_e = i

    # Initialize the aligned strings as the input strings up to the position of the high score.
    v_aligned, w_aligned = v[:i], w[:j]

    # Quick lambda function to insert indels.
    insert_indel = lambda word, i: word[:i] + '-' + word[i:]

    # Backtrack to start of the fitting alignment.
    while i*j != 0:
        if backtrack[i][j] == 0:
            i -= 1
            w_aligned = insert_indel(w_aligned, j)
        elif backtrack[i][j] == 1:
            j -= 1
            v_aligned = insert_indel(v_aligned, i)
        elif backtrack[i][j] == 2:
            i -= 1
            j -= 1

    v_s = i
    w_s = j
    # Cut off v at the ending point of the backtrack.
    v_aligned = v_aligned[i:]
    # if i == 0:
    #     print "{}\n{}\n{}".format(v, w, j)

    return AlignmentResult(v_s, v_e, w_s, w_e, v_aligned, w_aligned, max_score)

def is_valid_alignment(aln_result, error_rate, min_len):
    pass

def fitting_alignment2(v, w, k, m, e):

    aln = fitting_alignment(v, w[:m])
    ed = hamming(aln.v_aligned, aln.w_aligned)

    if ed > 1: return None

    M = [[0 for j in xrange(len(w)+1)] for i in xrange(len(v)+1)]

    for i in xrange(1,len(v)+1):
        M[i][0] = i
    for j in xrange(1,len(w)+1):
        M[0][j] = j

    for i in xrange(1,len(v)+1):
        for j in xrange(1,len(w)+1):
            if v[i-1] == w[j-1]:
                M[i][j] = M[i-1][j-1]
            else:
                M[i][j] = min(M[i-1][j]+1, M[i][j-1]+1, M[i-1][j-1]+1)

    i = aln.v_e
    j = m

    v_aligned = aln.v_aligned
    w_aligned = aln.w_aligned

    while (j+1 <= len(w) or i+1 <= len(v)):
        error_rates = [M[i+1][j]/float(j), M[i+1][j+1]/float(j+1), M[i][j+1]/float(j+1)]
        min_e = min(error_rates)
        if min_e > e:
            return AlignmentResult(aln.v_s, i, aln.w_s, j, v_aligned, w_aligned, 0)
        ind = error_rates.index(min_e)
        if ind == 0:
            i += 1
            v_aligned += v[i]
            w_aligned += '-'
        elif ind == 1:
            i += 1
            v_aligned += v[i]
            j += 1
            w_aligned += w[j]
        else:
            v_aligned += '-'
            j += 1
            w_aligned += w[j]

class TargetSequence(object):
    """docstring for TargetSequence"""
    def __init__(self, l_pos, seq):
        self.l_pos = l_pos
        self.seq = seq
        

class FiveEndSRead(object):
    """docstring for FiveEndSRead"""
    def __init__(self, id, l_pos, tail, seq):
        self.id = id
        self.l_pos = l_pos
        self.tail = tail
        self.seq = seq

    def __len__(self):
        return len(seq)

class FiveEndForwardSRead(FiveEndSRead):
    """docstring for FiveEndForwardSRead"""
    def __init__(self, id, l_pos, tail, seq):
        super(FiveEndForwardSRead, self).__init__(id, l_pos, tail, seq)

    def s_pos(self):
        return self.l_pos + self.tail

    def s_seq(self):
        return self.seq[:self.tail]

    def find_call(self, target_seq, error_rate):
        v = self.target_seq.seq
        w = self.read.s_seq()
        aln = fitting_alignment(v, w)
        if float(aln.edit_dist)/len(w) <= error_rate:
            return Call()

class NullCall(object):
    pass

class Call(object):
    """docstring for Call"""
    def __init__(self, bp1, bp2):
        self.bp1 = bp1
        self.bp2 = bp2
        
class ReadTargetSequence(object):
    """docstring for ReadTargetSequence"""
    def __init__(self, read, target_seq):
        self.read = read
        self.target_seq = target_seq


class AlignmentResult(object):
    """docstring for AlignmentResult"""
    def __init__(self, v_s, v_e, w_s, w_e, v_aligned, w_aligned, score, edit_dist):
        self.v_s = v_s
        self.v_e = v_s
        self.w_s = w_s
        self.w_e = w_e
        self.v_aligned = v_aligned
        self.w_aligned = w_aligned
        self.score = score
        self.edit_dist = edit_dist

    def error_rate(self):
        pass

    def alignment_length(self):
        pass
        
def fetch_target_seq(region):
    pass

def sp_region_forward(s_pos):
    pass

def sp_region_reverse(s_pos):
    pass

def fetch_m_pos_list(sp_region):
    pass

def focal_region(m_pos):
    pass

def merge_focal_regions(regions):
    pass

if __name__ == '__main__':

    # Read the input data.
    with open('data/chr15.132.dels.4x.svseq2.txt') as input_data, \
    open('output/chr15.132.dels.4x.svseq2.txt', 'w') as output_data:
        for line1, line2 in itertools.izip_longest(*[input_data]*2):
            # print(line1,line2)
            word1, word2 = line1.strip().upper(), line2.strip()

            # Get the fitting alignment.
            alignment = fitting_alignment(word1, word2)

            # Print and save the answer.
            # print '\n'.join(alignment)
            if int(alignment[0]) <= 2:
                output_data.write('\n'.join(itertools.imap(str, alignment))+'\n')
