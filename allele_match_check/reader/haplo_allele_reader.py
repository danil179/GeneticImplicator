# pylint: disable=C0301

'''
Contains haplo allele reader class
'''

from reader.haplo_reader import HaploReader
from genetic.haplo_sample import HaploSample


class HaploAlleleReader(HaploReader):
    '''
    Reader for haplotypes represented as strings of alleles
    '''

    def __init__(self, allele_names, name_seperator, allele_seperator):
        self._allele_names = allele_names
        self._order = dict((allele_name, i) for i, allele_name in enumerate(allele_names))
        self._allele_seperator = allele_seperator
        self._name_seperator = name_seperator

    def read_haplo_from_string(self, haplo_str):
        majorlis = [None] * len(self._allele_names)
        minorlis = [None] * len(self._allele_names)
        arr = haplo_str.split(self._allele_seperator)
        for allele in arr:
            allele_name = allele.split(self._name_seperator)[0].lstrip("'")
            major = int(allele.split(self._name_seperator)[1].split(':')[0])
            minor = int(''.join([i for i in allele.split(self._name_seperator)[1].split(':')[1] if i.isdigit()]))
            majorlis[self._order[allele_name]] = major
            minorlis[self._order[allele_name]] = [minor]
        return HaploSample(majorlis, minorlis, None, None)
