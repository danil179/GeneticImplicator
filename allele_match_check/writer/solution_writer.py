'''
Contains solution writer class
'''
import itertools
from genetic.parent_type import ParentType


class SolutionWriter:
    '''
    Writer for the solutions
    '''

    def __init__(self, filename, haplo_writer):
        self._filename = filename
        self._haplo_writer = haplo_writer
        self._init_file()

    def _init_file(self):
        self._csv_file = open(self._filename, 'w')
        self._csv_file.write('FAMCODE,F1,F2,M1,M2,CHILDS\n')
        self._csv_file.flush()

    def write_solution(self, solutions, family):
        """
        Write solutions to already opened csv file with child associations

        Args:
            solutions (list of MaxSolutionState): solutions to write
            family (Family): family

        """
        # many solution exist, iterating every solution
        for max_sol_state in solutions:
            # max_sol still has many solutions
            for max_sol in max_sol_state.get_state_list():
                self._write_line(max_sol, family, max_sol_state.get_free_haplos())

    def _write_line(self, max_sol, family, free_haplos):
        """Short summary.

        Args:
            max_sol (list of haplotypes): 4 haplotypes for parents of the solution
            family (Family): family
            free_haplos (dict[parent_type][hap_num] or booleans): True if the allele is free/empty

        """
        # write famcode
        stri = str(family.code)
        stri += ','
        # write parents (already ordered F1,F2,M1,M2)
        for i in range(4):
            stri += '\"'
            stri += self._haplo_writer.write_haplotype(max_sol[i])
            stri += '\"'
            stri += ','

        stri += '\"'

        # declare the print names for this order - introduces in max solution
        parents_names_index = ['F1', 'F2', 'M1', 'M2']
        ptype_index = {0: ParentType.FATHER, 1: ParentType.FATHER, 2: ParentType.MOTHER, 3: ParentType.MOTHER}
        haplo_index = {0: 0, 1: 1, 2: 0, 3: 1}
        # TODO: in the unordered case this is no longer true that inspecting one haplotype will give you solution, what we need is
        # to check child covering in all 4 possible coverings (F1,M1 - F1,M2 - F2,M1 - F2,M2), use match_unordered_haplos in human_sample
        # write childs, long nested for
        # e.g. 1=M1~F2:2=F1+M1~F1+M2
        for ch_num, child in enumerate(family.full_ordered_childs):
            if ch_num != 0:
                stri += ':'
            stri += str(family.child_indexes[ch_num])
            stri += '='
            # order corresponding to the index for each parent haplotype (order is F1,F2,M1,M2)
            # check all the matches of parents haplotypes with current child's haplotype, there are 4 options (F1,M1),(F1,M2),(F2,M1),(F2,M2)
            for orders in [[0, 2], [0, 3], [1, 2], [1, 3]]:
                if child.match_unordered_haplos_samp(max_sol[orders[0]], max_sol[orders[1]]):
                    stri += parents_names_index[orders[0]]
                    stri += '~'
                    stri += parents_names_index[orders[1]]
                    stri += '+'
            # removes last +
            stri = stri.rstrip('+')

        # DEBUG: not needed currently
        # stri += '\",\"'
        # save used parents string
        # for ptype in ParentType:
        #    for hapnum in range(2):
        #        if not free_haplos[ptype][hapnum]:
        #            stri += ptype.value + str(hapnum + 1) + "~"
        # remove last '~'
        # stri = stri.rstrip('~')
        stri += '\"\n'
        self._csv_file.write(stri)
        # REVIEW: consider buffering, depends on debug mode
        self._csv_file.flush()
