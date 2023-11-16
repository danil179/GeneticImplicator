'''
Contains allele types enum
'''

from enum import Enum


class AlleleType(Enum):
    '''
    Enum for the different allele parts (major/minor)
    '''
    MAJOR = 'Major'
    MINOR = 'Minor'
