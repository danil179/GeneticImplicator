'''
Contains haplo writer class.
'''

from writer.haplo_writer import HaploWriter
# REVIEW: consider merge this with haplo reader to avoid __init__ duplicated code


class HaploAlleleWriter(HaploWriter):
    '''
    Interface to write haplotypes
    '''

    def __init__(self, allele_names, name_seperator, allele_seperator):
        self._allele_names = allele_names
        self._order = dict((i, allele_name) for i, allele_name in enumerate(allele_names))
        self._allele_seperator = allele_seperator
        self._name_seperator = name_seperator

    def write_haplotype(self, haplo):
        """
        converts haplotype to string

        Args:
            haplo (HaploSample): haplotype to convert

        """
        stri = ''
        for i in range(haplo.allele_count):
            if i != 0:
                stri += self._allele_seperator
            stri += self._order[i]  # write allele name
            stri += self._name_seperator
            stri += str(haplo.majorli[i]).zfill(2)
            stri += ':'
            if len(haplo.minorli[i]) == 1:
                stri += str(haplo.minorli[i][0]).zfill(2)
            else:
                stri += str(haplo.minorli[i]).zfill(2)
        return stri
