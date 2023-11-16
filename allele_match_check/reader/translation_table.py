'''
Constains translation table class
'''


import logging


class TranslationTable:
    """
    Translate ambiguity string to list of alleles
    """

    ERROR_AMBIGUITY_EXISTS = 'Ambiguity already exists'

    def __init__(self):
        self._trans_table = dict()

    def resolve_minor(self, minor_str, major_allele):
        """
        Find the list of alleles that represents the ambiguities in the minor allele

        :param minor_str: minor's ambiguity string
        :type minor_str: string
        :param major_allele: major allele
        :type major_allele: int
        """
        if minor_str not in self._trans_table:  # ambiguity unresolved
            logging.info('missing ambiguity type(%s) - considered blank', minor_str)
            return ['']
        minor = self._trans_table[minor_str]
        if isinstance(minor, dict):  # different options for each major allele
            if major_allele not in minor:
                logging.info('inconsistent ambiguity type (%s-%d)', minor_str, major_allele)
                return None
            # in the serialization file there are duplicates in the arrays,
            # to remove them there is a set function and then sort to get only one order
            return sorted(set(minor[major_allele].tolist()))
        return sorted(set(minor.tolist()))

    def add_minor_ambiguity(self, ambiguity_str, alleles, major_allele=None):
        """
        Add ambiguity to the translation table, if already exists raise RuntimeError
        :param minor_str: ambiguity
        :type minor_str:
        :param alleles: list of the alleles for this ambiguity
        :type alleles: list/array
        :param major_allele:
        :type major_allele: int
        """

        if major_allele is None:
            if ambiguity_str in self._trans_table:
                raise RuntimeError(self.ERROR_AMBIGUITY_EXISTS)
            self._trans_table[ambiguity_str] = alleles
        else:
            if ambiguity_str not in self._trans_table:
                self._trans_table[ambiguity_str] = dict()
            if major_allele in self._trans_table[ambiguity_str]:
                raise RuntimeError(self.ERROR_AMBIGUITY_EXISTS)
            self._trans_table[ambiguity_str][major_allele] = alleles
