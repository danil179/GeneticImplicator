# pylint: disable=C0301

'''
Contains ambiguity reader class
'''

from abc import abstractmethod
import array
import collections
from reader.translation_table import TranslationTable


class AmbiguityReader:
    """Read a file that describes the types of ambiguities in data."""

    ERROR_WRONG_IMPLEMENT = 'interface AbiguityReader.GetAmbiguityAlleles implemented wrong'

    def __init__(self, data):
        """
        Create new ambiguity reader

        :param data: variable that holds the data (or a way to get the data)
        :type data:
        """
        self._data = data
        self.translation = TranslationTable()

    @abstractmethod
    def _get_next_entry(self):
        """
        Return a list of the alleles if there is no major ambiguity, or dictionary where the key is the major allele and the value is the list of the alleles, if the last entry reached returns None.

        :param entry_data: entry of one ambiguity (contains the data of the ambiguity)
        :type entry_data:
        """
        pass

    @abstractmethod
    def _get_ambiguity_alleles(self, entry_data):
        """
        Return a list of the alleles if there is no major ambiguity, or dictionary where the key is the major allele and the value is the list of the alleles.

        :param entry_data: entry of one ambiguity (contains the data of the ambiguity)
        :type entry_data:
        """
        pass

    @abstractmethod
    def _get_ambiguity_key(self, entry_data):
        """
        Return a key that represents the ambiguity.

        :param entry_data: entry of one ambiguity (contains the data of the ambiguity)
        :type entry_data:
        """
        pass

    @abstractmethod
    def _free_resources(self):
        """
        Free all the resources used to read ambiguities.
        """

    def read_ambiguities(self):
        """
        Reads the ambiguities and add that to the table.
        """
        # reads all ambiguities until the end of entries
        while True:
            entry = self._get_next_entry()
            if entry is None:
                self._free_resources()
                return
            key = self._get_ambiguity_key(entry)
            val = self._get_ambiguity_alleles(entry)
            if isinstance(val, (collections.Sequence, array.array)):
                self.translation.add_minor_ambiguity(key, val)
            elif isinstance(val, dict):
                for major in val.keys():
                    self.translation.add_minor_ambiguity(key, val[major], major)
            else:
                raise RuntimeError(self.ERROR_WRONG_IMPLEMENT)
