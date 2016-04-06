from abc import ABCMeta, abstractmethod
import re


def natural_key(s):
    return [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', s)]


class VariantTypes(object):
    """docstring for VariantType"""
    DEL = 0
    INS = 1
    INV = 2
    DUP = 3


class GenomePosition(object):
    """docstring for GenomePosition"""
    def __init__(self, ref_name, pos):
        self.ref_name = ref_name
        self.pos = pos

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.ref_name, self.pos) == (other.ref_name, other.pos)
        return NotImplemented

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            return (natural_key(self.ref_name), self.pos) < \
                (natural_key(other.ref_name), other.pos)
        return NotImplemented

    def __hash__(self):
        return hash((self.ref_name, self.pos))

    def __str__(self):
        return '{}:{}'.format(self.ref_name, self.pos)

    def genome_position_with_ci(self, slop):
        return GenomePositionWithCi(self, Interval(-slop, slop))


class Interval(object):
    """docstring for Interval"""
    def __init__(self, a, b):
        assert a <= b
        self.a = a
        self.b = b

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.a, self.b) == (other.a, other.b)
        return NotImplemented

    def __len__(self):
        return self.b - self.a

    def __str__(self):
        return '[{}, {}]'.format(self.a, self.b)

    def overlaps(self, other):
        return (self.a >= other.a and self.a < other.b) or \
            (other.a >= self.a and other.a < self.b)

    def reciprocal_overlaps(self, other, threshold):
        if not self.overlaps(other):
            return False

        o = min(self.b, other.b) - max(self.a, other.a)
        n = len(self)
        m = len(other)

        return float(o)/max(n, m) >= threshold


class GenomePositionWithCi(object):
    """docstring for GenomePositionWithCi"""
    def __init__(self, genome_pos, confidence_interval):
        self.genome_pos = genome_pos
        self.confidence_interval = confidence_interval

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.genome_pos, self.confidence_interval) == (other.genome_pos, other.confidence_interval)
        return NotImplemented

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            return self.genome_pos < other.genome_pos
        return NotImplemented

    def __str__(self):
        return '{}, {}'.format(self.genome_pos, self.confidence_interval)

    def __hash__(self):
        return hash(self.genome_pos)

    def to_genome_region(self):
        return GenomeRegion(
            GenomePosition(self.genome_pos.ref_name, self.genome_pos.pos + self.confidence_interval.a),
            GenomePosition(self.genome_pos.ref_name, self.genome_pos.pos + self.confidence_interval.b)
        )


class GenomeRegion(object):
    """docstring for GenomeRegion"""
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.start, self.end) == (other.start, other.end)
        return NotImplemented

    def overlaps(self, other):
        if self.start.ref_name != other.start.ref_name or self.end.ref_name != other.end.ref_name:
            return False
        return Interval(self.start.pos, self.end.pos).overlaps(Interval(other.start.pos, other.end.pos))

    def reciprocal_overlaps(self, other, threshold):
        if self.start.ref_name != other.start.ref_name or self.end.ref_name != other.end.ref_name:
            return False
        return Interval(self.start.pos, self.end.pos).reciprocal_overlaps(Interval(other.start.pos, other.end.pos), threshold)


class Variant(object):
    """docstring for Variant"""

    __metaclass__ = ABCMeta

    def __init__(self, name, pos1, pos2, micro_hom_seq=None, micro_ins_seq=None):
        self.name = name
        self.pos1 = pos1
        self.pos2 = pos2
        self.micro_hom_seq = micro_hom_seq
        self.micro_ins_seq = micro_ins_seq

    def __eq__(self, other):
        return (self.variant_type(), self.pos1, self.pos2) == \
            (other.variant_type(), other.pos1, other.pos2)

    def __lt__(self, other):
        return (self.pos1, self.pos2) < (other.pos1, other.pos2)

    def __hash__(self):
        return hash((self.variant_type(), self.pos1, self.pos2))

    def __str__(self):
        return '{}, {}, {}, L:{}'.format(self.pos1, self.pos2, self.name, len(self))

    def __len__(self):
        return self.pos2.genome_pos.pos - self.pos1.genome_pos.pos - 1

    def _genome_region(self):
        return GenomeRegion(self.pos1.genome_pos, self.pos2.genome_pos)

    def matched_with_both_overlaps(self, other):
        if self.variant_type() != other.variant_type():
            return False
        return self.pos1.to_genome_region().overlaps(other.pos1.to_genome_region()) and \
            self.pos2.to_genome_region().overlaps(other.pos2.to_genome_region())

    def matched_with_reciprocal_overlap(self, other, threshold):
        if self.variant_type() != other.variant_type():
            return False
        return self._genome_region().reciprocal_overlaps(other._genome_region(), threshold)

    @abstractmethod
    def variant_type(self):
        pass


class Deletion(Variant):
    """docstring for Deletion"""

    def variant_type(self):
        return VariantTypes.DEL


class VariantFile(object):
    """docstring for VariantFile"""

    __metaclass__ = ABCMeta

    def __init__(self, filename):
        self.filename = filename

    def __enter__(self):
        self._file_obj = self._open_file()
        return self

    def __exit__(self, *args):
        self._close_file()

    def __iter__(self):
        for line in self._file_obj:
            if self.is_valid_line(line):
                data = self.line_to_data(line)
                variant_type = self.get_variant_type(data)
                if variant_type == VariantTypes.DEL:
                    yield self.data_to_deletion(data)
                if variant_type == VariantTypes.INS:
                    yield self.data_to_insertion(data)
                if variant_type == VariantTypes.INV:
                    yield self.data_to_inversion(data)
                if variant_type == VariantTypes.DUP:
                    yield self.data_to_duplication(data)

    def _open_file(self):
        return open(self.filename, 'r')

    def _close_file(self):
        self._file_obj.close()

    @abstractmethod
    def get_variant_type(self, data):
        pass

    def is_valid_line(self, line):
        return True

    def line_to_data(self, line):
        return line.split()

    @abstractmethod
    def data_to_deletion(self, data):
        pass

    @abstractmethod
    def data_to_insertion(self, data):
        pass

    @abstractmethod
    def data_to_inversion(self, data):
        pass

    @abstractmethod
    def data_to_duplication(self, data):
        pass


class SvsimBedpeFile(VariantFile):
    """docstring for SvsimBedpeFile"""

    def __init__(self, filename, slop):
        super(SvsimBedpeFile, self).__init__(filename)
        self.slop = slop

    def get_variant_type(self, data):
        variant_types = ['DEL', 'INS', 'INV', 'DUP']
        pattern = '(?P<variant_type>{})'.format('|'.join(variant_types))
        m = re.search(pattern, data[6])
        return variant_types.index(m.group('variant_type'))

    def data_to_deletion(self, data):
        return Deletion(
            data[6],
            GenomePosition(data[0], int(data[2])).genome_position_with_ci(self.slop),
            GenomePosition(data[3], int(data[5])).genome_position_with_ci(self.slop)
        )

    def data_to_insertion(self, data):
        pass

    def data_to_inversion(self, data):
        self._file_obj.next()
        return Deletion(
            data[6],
            GenomePosition(data[0], int(data[2])),
            GenomePosition(data[3], int(data[5]))
        )

    def data_to_duplication(self, data):
        pass


class SpritesBedpeFile(VariantFile):
    """docstring for SpritesBedpeFile"""

    def get_variant_type(self, data):
        variant_types = ['DEL', 'INS', 'INV', 'DUP']
        pattern = '(?P<variant_type>{})'.format('|'.join(variant_types))
        m = re.search(pattern, data[6])
        return variant_types.index(m.group('variant_type'))

    def data_to_deletion(self, data):
        delta = int(data[2]) - int(data[1]) - 1
        if data[6].endswith('5F'):
            return Deletion(
                data[6],
                GenomePositionWithCi(GenomePosition(data[0], int(data[2])), Interval(-delta, 0)),
                GenomePositionWithCi(GenomePosition(data[3], int(data[5])), Interval(-delta, 0))
            )
        if data[6].endswith('5R'):
            return Deletion(
                data[6],
                GenomePositionWithCi(GenomePosition(data[0], int(data[1])+1), Interval(0, delta)),
                GenomePositionWithCi(GenomePosition(data[3], int(data[4])+1), Interval(0, delta))
            )
        return Deletion(
            data[6],
            GenomePositionWithCi(GenomePosition(data[0], int(data[2])), Interval(-delta, 0)),
            GenomePositionWithCi(GenomePosition(data[3], int(data[5])), Interval(-delta, 0))
        )

    def data_to_insertion(self, data):
        pass

    def data_to_inversion(self, data):
        pass

    def data_to_duplication(self, data):
        pass
