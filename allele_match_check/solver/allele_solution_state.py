# pylint: disable=C0301
'''
Contains allele solution state class
'''


import copy
from genetic.parent_type import ParentType
from genetic.allele_type import AlleleType


class AlleleSolutionState:
    """
    represents a state in generating the solutions for the single allele problem
    """

    def __init__(self):
        self._allele_state = dict()

        self._order_index = 0

        for ptype in ParentType:
            self._allele_state[ptype] = dict()
            for atype in AlleleType:
                self._allele_state[ptype][atype] = [[], []]

    @classmethod
    def create_from_haplo(cls, haplo, allele_ind):
        """
        Initiate state from haplotype

        Args:
            haplo (HaploSample): the haplotype
            allele_ind (int): index of the allele to start the state from

        """
        state = cls()
        state.add_allele(haplo.majorli[allele_ind], haplo.minorli[allele_ind], haplo.haptype, haplo.hapnum)
        return state

    def add_haplotype(self, haplo):
        if haplo.haptype is None or haplo.hapnum is None:
            raise Exception("haplotype type/num cannot be nulled")
        for allele_ind in range(len(haplo.majorli)):
            self.add_allele(haplo.majorli[allele_ind], haplo.minorli[allele_ind], haplo.haptype, haplo.hapnum)

    def __add__(self, other):
        # this adding alleles assuming internal order, without order one should consider all 4 possibilities of hap_num and ptype
        res = AlleleSolutionState()
        for ptype in ParentType:
            for atype in AlleleType:
                for hap_num in range(2):
                    res._allele_state[ptype][atype][hap_num] = self._allele_state[ptype][atype][hap_num] + other._allele_state[ptype][atype][hap_num]  # pylint: disable=W0212
        return res

    def add_ignoring_order(self, other):
        # init result
        res = [None] * 4

        # 4 options
        for i, order in enumerate([[False, False], [False, True], [True, False], [True, True]]):
            other1 = copy.copy(other)
            if order[0]:
                other1.swap_parent_haplotype(ParentType.MOTHER)
            if order[1]:
                other1.swap_parent_haplotype(ParentType.FATHER)
            res[i] = self + other1
        return res

    def swap_parent_haplotype(self, ptype):
        for atype in AlleleType:
            # swap hap_num
            temp = self._allele_state[ptype][atype][1]
            self._allele_state[ptype][atype][1] = self._allele_state[ptype][atype][0]
            self._allele_state[ptype][atype][0] = temp

    def check_permutation_state(self, permutation_state):
        """
        Returns True if permutation_state is a permutation of current state that covered by current state.
        :param permutation_state: the state to check
        :type permutation_state: AlleleSolutionState
        """
        for porder in [[ParentType.FATHER, ParentType.MOTHER], [ParentType.MOTHER, ParentType.FATHER]]:
            for haporder in [[0, 0], [0, 1], [1, 0], [1, 1]]:
                # haporder[0] = if father haplotypes swapped
                # haporder[1] = if mother haplotypes swapped

                is_permutation = True
                # check if the swap cause covering
                for ptype in ParentType:
                    if ptype == ParentType.FATHER:
                        swap_index = 0
                    elif ptype == ParentType.MOTHER:
                        swap_index = 1
                    # long if to check the 4 possible swaps in current parent, possible to implement with 2 for loops (one for allele type, and one for haplotype number)
                    if not (self._compare_allele_list(self.get_all_alleles(ptype, AlleleType.MAJOR, 0), permutation_state.get_all_alleles(porder[swap_index], AlleleType.MAJOR, 0 + haporder[swap_index]))
                            and self._compare_allele_list(self.get_all_alleles(ptype, AlleleType.MAJOR, 1), permutation_state.get_all_alleles(porder[swap_index], AlleleType.MAJOR, 1 - haporder[swap_index]))
                            and self._compare_allele_list(self.get_all_alleles(ptype, AlleleType.MINOR, 0), permutation_state.get_all_alleles(porder[swap_index], AlleleType.MINOR, 0 + haporder[swap_index]))
                            and self._compare_allele_list(self.get_all_alleles(ptype, AlleleType.MINOR, 1), permutation_state.get_all_alleles(porder[swap_index], AlleleType.MINOR, 1 - haporder[swap_index]))):
                        is_permutation = False
                        break
                if is_permutation:
                    return True
        return False

    @staticmethod
    def _compare_allele_list(covering_allist, allist):
        """
        Returns True if covering_allist covers allist s.t allist is a subset of covering_allist.
        :param covering_allist: the covering allele list to check
        :type covering_allist: list
        :param allist: covered allele list
        :type allist: the covered allele list to check
        """
        for i, allele in enumerate(covering_allist):
            # for majors al[i] is an int, for minor it's a list of ints
            # minor
            if isinstance(allele, list):
                if allele != [''] and not set(allist[i]).issubset(set(allele)):
                    return False
            # major
            else:
                if allele not in (allist[i], ''):
                    return False
        return True

    def compare_state(self, compare_state):
        """
        return the state that is covered (subset of the other state), if no state cover the other return None.
        if the states are equal then return compare_state
        :param compare_state: state to compare with
        :type compare_state: AlleleSolutionState
        """
        compare_state_covered = True
        current_state_covered = True
        for ptype in ParentType:
            for atype in AlleleType:
                for hapnum in range(2):
                    # instead calling to _CompareAlleleList twice (then it will take 2 loops), using compare in parallel
                    for al_ind, al_list in enumerate(self.get_all_alleles(ptype, atype, hapnum)):
                        # two options: major allele, minor allele
                        # minor allele
                        compare_allele = compare_state.get_allele(ptype, atype, hapnum, [al_ind])[0]
                        if atype == AlleleType.MINOR:
                            # al_list = [5,6,7,10]
                            # match only if the new solution covered by existing solutions
                            compare_state_covered = compare_state_covered and (al_list == [''] or set(compare_allele).issubset(al_list))
                            current_state_covered = current_state_covered and (compare_allele == [''] or set(al_list).issubset(compare_allele))
                        # major allele
                        else:
                            # checking if each allele is subset of the other and assign the result to the coresponding variable
                            compare_state_covered = compare_state_covered and al_list in ('', compare_allele)
                            current_state_covered = current_state_covered and compare_allele in ('', al_list)
                        if not compare_state_covered and not current_state_covered:
                            return None
        if compare_state_covered:
            return compare_state
        if current_state_covered:
            return self
        return None

    def _add_major_allele(self, allele, parent_type, hap_num):
        """
        Adding major allele to the state
        :param allele:
        :type allele:
        :param parent_type:
        :type parent_type:
        :param hap_num:
        :type hap_num:
        """

        self._allele_state[parent_type][AlleleType.MAJOR][hap_num].append(allele)

    def _add_minor_allele(self, alleles, parent_type, hap_num):
        """
        Adding minor allele to the state
        :param alleles:
        :type alleles:
        :param parent_type:
        :type parent_type:ParentType
        :param hap_num:
        :type hap_num:
        """
        self._allele_state[parent_type][AlleleType.MINOR][hap_num].append(alleles)

    def update_major_allele(self, allele, parent_type, hap_num, allele_num=0):
        """
        Update major part of the allele, if not exist raise exception

        Args:
            allele (int): the new major part
            parent_type (ParentType): parent type
            hap_num (type): haplotype number
            allele_num (type): index of allele to update

        """
        self._allele_state[parent_type][AlleleType.MAJOR][hap_num][allele_num] = allele
        return self

    def update_minor_allele(self, alleles, parent_type, hap_num, allele_num=0):
        """
        Update minor part of the allele, if not exist raise exception

        Args:
            alleles (list of int): new ambiguities in the allele
            parent_type (ParentType): parent type
            hap_num (int): haplotype number
            allele_num (int): index of allele to update

        """
        self._allele_state[parent_type][AlleleType.MINOR][hap_num][allele_num] = alleles
        return self

    def update_allele(self, major, minorlist, parent_type, hap_num, allele_num=0):
        """
        Update allele, if not exist raise exception

        Args:
            major (int): the new major part
            minorlist (list of int): new ambiguities in the allele
            parent_type (ParentType): parent type
            hap_num (int): haplotype number
            allele_num (int): index of allele to update

        """
        self.update_minor_allele(minorlist, parent_type, hap_num, allele_num)
        self.update_major_allele(major, parent_type, hap_num, allele_num)
        return self

    def add_allele(self, major, minorlist, parent_type, hap_num):
        """
        Add allele to the state, appending it the exisiting alleles if exist.

        Args:
            major (int): major part of the allele
            minorlist (list of int): minor part of the allele
            parent_type (ParentType): type of the parent
            hap_num (int): haplotype number (0/1)

        """
        self._add_minor_allele(minorlist, parent_type, hap_num)
        self._add_major_allele(major, parent_type, hap_num)
        return self

    def increase_order_index(self):
        """
        Icreasing the index of the order (that applies when solving for single allele, and keeping the order in the iteration)
        """
        self._order_index += 1
        # allow chaining
        return self

    def get_order_index(self):
        """
        Returns the current order index (see @increase_order_index)
        """
        return self._order_index

    def get_allele(self, ptype, atype, hapnum, order_allele_ind_list=None):
        """
        Returns a list of the specific alleles that requested
        :param ptype:
        :type ptype: ParentType
        :param atype:
        :type atype: AlleleType
        :param hapnum: haplotype number
        :type hapnum: int
        :param order_allele_ind: allele index to get according to the internal order
        :type order_allele_ind: int
        """
        if order_allele_ind_list is None:
            order_allele_ind_list = [0]
        # if no index exist for that allele return empty list
        if any(len(self._allele_state[ptype][atype][hapnum]) <= allele_ind for allele_ind in order_allele_ind_list):
            return []
        return [self._allele_state[ptype][atype][hapnum][order_allele_ind] for order_allele_ind in order_allele_ind_list]

    def get_all_alleles(self, ptype, atype, hapnum):
        """
        Returns all alleles (as ordered list) at this specific position
        :param ptype:
        :type ptype: ParentType
        :param atype:
        :type atype: AlleleType
        :param hapnum: haplotype number
        :type hapnum: int
        """
        return self._allele_state[ptype][atype][hapnum]

    def __copy__(self):
        new_state = AlleleSolutionState()
        new_state._allele_state = self.__get_state()  # pylint: disable=W0212
        new_state._order_index = self._order_index  # pylint: disable=W0212
        return new_state

    def __repr__(self):
        return str(self)

    def __str__(self):
        return str(self._allele_state)

    def __get_state(self):
        """
        Returns representation of the state as hash table (internal variable)
        """
        return copy.deepcopy(self._allele_state)

    def __hash__(self):
        return hash(tuple([tuple([tuple(str(allele)) for allele in self._allele_state[ptype][atype][hapnum]]) for ptype in ParentType for atype in AlleleType for hapnum in range(2)]))

    def __eq__(self, other):
        ''' minor list is a list of lists, checking equality while ignoring order, but because all the orders of the same list should be same then only check equality
        '''
        for ptype in ParentType:
            for atype in AlleleType:
                if self._allele_state[ptype][atype] != other._allele_state[ptype][atype]:  # pylint: disable=W0212
                    return False
        return True
