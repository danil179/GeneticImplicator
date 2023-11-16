'''
Module for parent types
'''
from enum import Enum


class ParentType(Enum):
    '''
    Enum for the different parents (mother/father)
    '''
    # Don't change the order, important for some iterations
    FATHER = 'F'
    MOTHER = 'M'
