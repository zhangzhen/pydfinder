import pysam

samfile = pysam.AlignmentFile("ex1.bam", "rb")

iter = samfile.fetch("seq1", 10, 20)
for x in iter:
    print (str(x))

class InsertSize(object):
	"""docstring for InsertSize"""
	def __init__(self, mean, sd, fold):
		self.mean = mean
		self.sd = sd
		self.fold = fold

	def len_in_normal_range(self, map_dist):
		return map_dist >= self.mean - self.fold * self.sd and map_dist <= self.mean + self.fold * self.sd

class SoftClippingRead(object):
	def __init__(self, id, gpos, map_dist):
		self.id = id
		self.gpos = gpos
		self.map_dist = map_dist

	def len_consistent(self, del, isize):
		if self.map_dist <= len(del):
			return False
		return isize.len_in_normal_range(self.map_dist - len(del))

class ThreePrimeEndForward(SoftClippingRead):
	"""docstring for ThreePrimeEndForward"""
	def __init__(self, gpos, map_dist):
		super(ThreePrimeEndForward, self).__init__(gpos, map_dist)

	def close_to(self, del, slop):
		return make_gregion(del.start(), slop).in_region(self.gpos)

class ThreePrimeEndReverse(SoftClippingRead):
	"""docstring for ThreePrimeEndReverse"""
	def __init__(self, gpos, map_dist):
		super(ThreePrimeEndReverse, self).__init__(gpos, map_dist)

	def close_to(self, del, slop):
		return make_gregion(del.end(), slop).in_region(self.gpos)
		

class GenomePosition(object):
	"""docstring for GenomePosition"""
	def __init__(self, chrom, pos):
		self.chrom = chrom
		self.pos = pos
		
def check_support(sc_read, del, slop, isize):
	return sc_read.len_consistent(del, isize) and sc_read.close_to(del, slop)
		
class BreakPoint(object):
	"""docstring for BreakPoint"""
	def __init__(self, left, right):
		self.gregion = GenomeRegion(left, right)

	def left(self):
		return self.gregion.start_gpos

	def right(self):
		return self.gregion.end_gpos
		
	def __len__(self):
		return len(gregion)

class Deletion(object):
	"""docstring for Deletion"""
	def __init__(self, id, breakpoint):
		self.id = id
		self.breakpoint = breakpoint

	def __len__(self):
		return len(self.breakpoint)

	def start(self):
		return self.breakpoint.left()

	def end(self):
		return self.breakpoint.right()

		
def make_gregion(gpos, slop):
	return GenomeRegion(GenomePosition(gpos.chrom, gpos.pos - slop),
		GenomePosition(gpos.chrom, gpos.pos + slop))

class GenomeRegion(object):
	
	def __init__(self, start_gpos, end_gpos):
		assert start_gpos.chrom == end_gpos.chrom and start_gpos.pos < end_gpos.pos
		self.start_gpos = start_gpos
		self.end_gpos = end_gpos

	def __len__(self):
		return self.end_pos - self.start_pos + 1

	def chromosome(self):
		return self.start_gpos.chrom

	def start_pos(self):
		return self.start_gpos.pos

	def end_pos(self):
		return self.end_gpos.pos

	def in_region(self, gpos):
		return gpos.chrom == self.chromosome() and gpos.pos >= self.start() and gpos.pos <= self.end()

class OverlapMatch(object):
	"""docstring for OverlapMatch"""
	def __init__(self, a, b, size, score, num_mismatches):
		self.a = a
		self.b = b
		self.size = size
		self.score = score
		self.num_mismatches = num_mismatches

def is_soft_clipping_read(read, min_qual):
	return read.is_proper_pair and read.mapping_quality >= min_qual		

