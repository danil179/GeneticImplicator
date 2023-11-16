'''
Contains haplo reader class.
'''

from abc import ABCMeta, abstractmethod


class HaploReader(metaclass=ABCMeta):
    '''
    Interface to read haplotypes
    '''
    @abstractmethod
    def read_haplo_from_string(self, haplo_str):
        '''
        Converts string to HaploSample
        :param haplo_str: string of the haplotype to convert
        :type haplo_str: string
        '''
        pass

    def read_haplotypes(self, filename):
        """
        Read each line as haplotype and return a list of all haplotypes.

        Args:
            filename (string): path to the file including name

        """
        filehaplo = open(filename, 'r')
        haplotypes = []
        for line in filehaplo.read().splitlines():
            haplotypes.append(self.read_haplo_from_string(line))
        return haplotypes
