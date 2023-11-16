'''
Contains data reader class.
'''

from abc import ABCMeta, abstractmethod


class DataReader(metaclass=ABCMeta):
    """Read the genetic data of the families"""

    def __init__(self, data, read_nulled_childs=False, remove_redundent=False):
        """
        Creates new reader
        :param data: variable that holds the information to get the data
        :type data:
        :param read_nulled_childs: if True then also read nulled childs (default False)
        :type read_nulled_childs: boolean
        :param remove_redundent: if True then removes redundant childs (default False)
        :type remove_redundent: boolean
        """
        self.data = data
        self.read_nulled_childs = read_nulled_childs
        self.remove_redundent = remove_redundent

    @abstractmethod
    def read_next_family(self):
        """Read the next family and returns Family"""
        pass

    def read_families(self):
        """Read all the families (ordered) from the current position and returns list"""
        families = []
        fam = self.read_next_family()
        while fam is not None:
            families.append(fam)
            fam = self.read_next_family()
        return families
