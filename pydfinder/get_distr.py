'''Calculate a distribution of insert sizes of read pairs
that span the clipping position of a soft-clipping read

Usage:
  get_distr.py -l <int> --is <mean_sd> <file> <position>
  get_distr.py -h | --help
  get_distr.py --version

Arguments:
  <file>                The BAM file
  <position>            A chromosome position MUST be specified as RNAME:POS

Options:
  -h --help                         Show this screen
  --version                         Show version
  -l <int>, --read-length <int>     Specify read length
  --is <mean_sd>                    Specify mean and sd as MEAN,SD
'''

from docopt import docopt
import pysam

__version__ = "0.1.0"
__author__ = "Zhen Zhang"
__license__ = "MIT"


class ChromosomePosition(object):
    """docstring for ChromosomePosition"""
    def __init__(self, chr_name, val):
        self.chr_name = chr_name
        self.val = val


class ChromosomeRegion(object):
    """docstring for ChromosomeRegion"""
    def __init__(self, chr_name, start, end):
        self.chr_name = chr_name
        self.start = start
        self.end = end


class LibraryParams(object):
    """docstring for LibraryParams"""
    def __init__(self, read_length, mean_insert_size, sd_insert_size):
        self.read_length = read_length
        self.mean_insert_size = mean_insert_size
        self.sd_insert_size = sd_insert_size

    def distance1(self):
        return self.read_length + self.mean_insert_size + 3*self.sd_insert_size

    def distance2(self):
        return self.read_length + self.mean_insert_size + 3*self.sd_insert_size

    def distance3(self):
        return 2*(self.mean_insert_size + 3*self.sd_insert_size)

def get_chr_region1(clipping_pos, lib_params):
    return ChromosomeRegion(clipping_pos.chr_name, clipping_pos.val, \
        clipping_pos.val + lib_params.distance1())

def get_chr_region2(clipping_pos, lib_params):
    return ChromosomeRegion(clipping_pos.chr_name, clipping_pos.val, \
        clipping_pos.val + lib_params.distance2())

def get_chr_region3(clipping_pos, lib_params):
    return ChromosomeRegion(clipping_pos.chr_name, clipping_pos.val, \
        clipping_pos.val + lib_params.distance3())

def is_both_mapped(read):
    return read.is_paired and not read.is_unmapped and not read.mate_is_unmapped \
        and read.reference_id == read.next_reference_id

def is_valid(read, clipping_pos):
    return is_both_mapped(read) and read.is_reverse and not read.mate_is_reverse \
        and read.reference_start >= clipping_pos \
        and read.next_reference_start + read.query_length <= clipping_pos

def get_insert_sizes(chr_region, samfile):
    return [abs(read.template_length) for read in samfile.fetch(chr_region.chr_name, chr_region.start-1, chr_region.end-1) \
        if is_valid(read, chr_region.start)]

def main():
    args = docopt(__doc__, version=__version__)
    samfile = pysam.AlignmentFile(args['<file>'], 'rb')
    chr_name, val = args['<position>'].split(':')
    clipping_pos = ChromosomePosition(chr_name, int(val))
    mean, sd = args['--is'].split(',')
    lib_params = LibraryParams(int(args['--read-length']), int(mean), int(sd))
    print "Normal distance:"
    print get_insert_sizes(get_chr_region1(clipping_pos, lib_params), samfile)
    print "Short distance:"
    print get_insert_sizes(get_chr_region2(clipping_pos, lib_params), samfile)
    print "Long distance:"
    print get_insert_sizes(get_chr_region3(clipping_pos, lib_params), samfile)


if __name__ == '__main__':
    main()
