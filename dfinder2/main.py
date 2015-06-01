import pysam

samfile = pysam.AlignmentFile("ex1.bam", "rb")

iter = samfile.fetch("seq1", 10, 20)
for x in iter:
    print (str(x))

class GenomeRegion(object):
	
	def __init__(self, chrom, start, end):
		self.chrom = chrom
		self.start = start
		self.end = end

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

