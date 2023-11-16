'''
Contains haplo writer class.
'''

from abc import ABCMeta, abstractmethod


class HaploWriter(metaclass=ABCMeta):
    '''
    Interface to write haplotypes
    '''

    @abstractmethod
    def write_haplotype(self, haplo):
        """
        converts haplotype to string

        Args:
            haplo (HaploSample): haplotype to convert

        """
        pass
