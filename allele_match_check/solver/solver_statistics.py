# pylint: disable=C0301
'''
Contains solver statistics class
'''


class SolverStatistics():
    '''
    Statistics for the FamilySolver class
    '''

    def __init__(self):
        # TODO: initialize all the statistics
        self.solved = 0
        self.nulled_in_single_allele = 0
        self.nulled_in_merge = 0
        self.nulled_in_max = 0
        self.nulled_in_haplo_check = 0

    @property
    def nulled(self):
        """
        Returns the number of nulled solutions
        """
        return self.nulled_in_single_allele + self.nulled_in_merge + self.nulled_in_max + self.nulled_in_haplo_check
