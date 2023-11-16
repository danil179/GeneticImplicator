# pylint: disable=C0301

'''
Contains max solution state class.
'''

import itertools
from genetic.parent_type import ParentType
from genetic.allele_type import AlleleType
from genetic.haplo_sample import HaploSample
from solver.allele_solution_state import AlleleSolutionState
import constants


class MaxSolutionState():
    """
    Final solution for the matching problem (maximum probability haplotype from haplotype list that consistent)
    """

    def __init__(self, allele_state):
        self._haplos = dict()
        for ptype in ParentType:
            self._haplos[ptype] = [[], []]

        self._max_freq = dict()
        for ptype in ParentType:
            # 0 probabilty for initialization
            self._max_freq[ptype] = [0, 0]

        self._haplo_free = dict()
        for ptype in ParentType:
            # which haplotype is free
            self._haplo_free[ptype] = [False, False]

        self._max_solution = dict()
        for ptype in ParentType:
            self._max_solution[ptype] = [[], []]

        self._load_haplotypes_from_allele_state(allele_state)

    def get_free_haplos(self):
        """
        Return dictionary dict, if dict[ptype][hapnum] is True then the haplotype is free
        and not constrained at all
        """
        return self._haplo_free

    def check_permutation_state(self, permutation_state):
        state1 = AlleleSolutionState()
        state2 = AlleleSolutionState()
        for ptype in ParentType:
            for hapnum in range(2):
                self._max_solution[ptype][hapnum][0].haptype = ptype
                self._max_solution[ptype][hapnum][0].hapnum = hapnum
                permutation_state._max_solution[ptype][hapnum][0].haptype = ptype
                permutation_state._max_solution[ptype][hapnum][0].hapnum = hapnum
                state1.add_haplotype(self._max_solution[ptype][hapnum][0])
                state2.add_haplotype(permutation_state._max_solution[ptype][hapnum][0])
        # now state1 and state2 are parent representation of this MaxSolutionState

        if state1.check_permutation_state(state2):
            return True
        return False

    def _load_haplotypes_from_allele_state(self, state):
        for ptype in ParentType:
            for hap_num in range(2):
                major = state.get_all_alleles(ptype, AlleleType.MAJOR, hap_num)
                minor = state.get_all_alleles(ptype, AlleleType.MINOR, hap_num)
                sample = HaploSample(major, minor, ptype, hap_num)
                self._haplos[ptype][hap_num] = sample
                self._haplo_free[ptype][hap_num] = sample.is_empty()
                if self._haplo_free[ptype][hap_num]:
                    self._max_freq[ptype][hap_num] = constants.TOTAL_FREQ
                    sample.haptype = None
                    sample.hapnum = None
                    self._max_solution[ptype][hap_num] = [sample]

    def _add_solution(self, haplo, ptype, hap_num):
        self._max_solution[ptype][hap_num].append(haplo)

    def _test_solution(self, haplo, ptype, hap_num):
        # if this haplo is free
        if self._haplo_free[ptype][hap_num]:
            return
        # not covering - not a solution
        if not self._haplos[ptype][hap_num].cover(haplo.haplotype):
            return
        # probability higher - remove current solution and update to the new one
        if haplo.frequency > self._max_freq[ptype][hap_num]:
            self._max_solution[ptype][hap_num] = [haplo.haplotype]
            self._max_freq[ptype][hap_num] = haplo.frequency
        # probability is the same - add the solution to the existing ones
        elif haplo.frequency == self._max_freq[ptype][hap_num]:
            self._add_solution(haplo.haplotype, ptype, hap_num)

    def process_haplo(self, haplo):
        """
        Check if the haplotype is better solution from previous ones.

        Args:
            haplo (HaploSample): the haplotype

        """
        for ptype in ParentType:
            for hap_num in range(2):
                self._test_solution(haplo, ptype, hap_num)

    def get_max_frequency(self):
        """
        Returns the current frequency (maximum among the already checked haplotypes)

        """
        max_freq = 1
        for ptype in ParentType:
            for hap_num in range(2):
                max_freq *= self._max_freq[ptype][hap_num]
        return max_freq

    def get_solution_haplotype(self, ptype, hap_num):
        """
        Returns current solution haplotypes as list of HaploSample (without ambiguities)

        Args:
            ptype (ParentType): parent type
            hap_num (int): haplotype number

        """
        return self._max_solution[ptype][hap_num]

    def get_state_list(self):
        """
        Returns list with all solutions ordered F1,F2,M1,M2

        """
        sol_list = [self._max_solution[ptype][hapnum] for ptype in ParentType for hapnum in range(2)]
        max_solution_list = list(itertools.product(*sol_list))
        return max_solution_list

    def __str__(self):
        stri = '\n'
        for ptype in ParentType:
            for hap_num in range(2):
                stri += ptype.value + str(hap_num + 1) + ': ' + str(self._max_solution[ptype][hap_num]) + '\n'
        return stri

    def __repr__(self):
        return str(self)
