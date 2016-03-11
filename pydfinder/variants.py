from abc import ABCMeta, abstractmethod
import re


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


class ConfidenceInterval(object):
    """docstring for ConfidenceInterval"""
    def __init__(self, a, b):
        self.a = a
        self.b = b


class ImpreciseGenomePosition(object):
    """docstring for ImpreciseGenomePosition"""
    def __init__(self, genome_pos, ci):
        self.genome_pos = genome_pos
        self.ci = ci


class GenomeRegion(object):
    """docstring for GenomeRegion"""
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.start, self.end) == (other.start, other.end)
        return NotImplemented


class Variant(object):
    """docstring for Variant"""

    __metaclass__ = ABCMeta

    def __init__(self, name, region, micro_hom_seq=None, micro_ins_seq=None):
        self.name = name
        self.region = region
        self.micro_hom_seq = micro_hom_seq
        self.micro_ins_seq = micro_ins_seq

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.name, self.region, self.micro_hom_seq, self.micro_ins_seq) == \
                (other.name, other.region, other.micro_hom_seq, other.micro_ins_seq)
        return NotImplemented

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

    def get_variant_type(self, data):
        variant_types = ['DEL', 'INS', 'INV', 'DUP']
        pattern = "r'(?P<variant_type>{})'".format('|'.join(variant_types))
        m = re.search(pattern, data[6])
        return variant_types.index(m.group('variant_type'))

    def data_to_deletion(self, data):
        return Deletion(data[6], GenomeRegion(GenomePosition(data[0], int(data[2])),
            GenomePosition(data[3], int(data[5]))), None, None)

    def data_to_insertion(self, data):
        pass

    def data_to_inversion(self, data):
        self._file_obj.next()
        return Deletion(data[6], GenomeRegion(GenomePosition(data[0], int(data[2])),
            GenomePosition(data[3], int(data[5]))), None, None)

    def data_to_duplication(self, data):
        pass


def get_variants(cls, filename):
    with cls(filename) as f:
        for variant in f:
            yield variant
