'''
Contains haplotype with probability class
'''


class HaploFreqSample():
    '''
    Describing haplotype sample with probability
    '''

    def __init__(self, haplo, freq):
        self._haplo = haplo
        self._freq = freq

    @property
    def frequency(self):
        """
        The frequency of the haplotype (higher frequency means higher probability)
        """
        return self._freq

    @property
    def haplotype(self):
        """
        The haplotype sample
        """
        return self._haplo
