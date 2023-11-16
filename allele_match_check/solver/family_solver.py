# pylint: disable=C0301

'''
Contains family solver class
'''

import copy
import logging
from _collections import deque
from genetic.parent_type import ParentType
from genetic.allele_type import AlleleType
from solver.allele_solution_state import AlleleSolutionState
from solver.solver_statistics import SolverStatistics
from solver.max_solution_state import MaxSolutionState


# TODO: consider splitting this class into general interface then abstract class of iterative solver (single solve, merge, and maximum) and 3 classes each handles different part of the solution


class FamilySolver:
    '''
    Family solver, solving the matching problem for families
    '''

    def __init__(self):
        self._statistics = SolverStatistics()

    def solve(self, family, haplotypes):
        """
        solves the general problem of finding the most likely parents
        and return the solution as FamilySolution

        :param family: family to solve
        :type family: Family
        :param haplotypes: all the haplotypes that have frequency and valid to the solutions
        :type haplotypes: HaploFreqSample
        """

        # FIXME: remove this or change (e.g. count different alleles num foreach allele)
        # if not self._fast_haplotype_check(family.childs):
        #    logging.debug('no solution found in first haplotype check')
        #    self._statistics.nulled_in_haplo_check += 1
        #    return None

        single_allele_sols = self._solve_match_in_alleles(family)
        if single_allele_sols is None:
            logging.debug('no solution in one of the alleles - continue to the next family')
            self._statistics.nulled_in_single_allele += 1
            return None

        # logging num of solutions in single allele
        maxnum = 0
        for i, single_sol in enumerate(single_allele_sols):
            newmax = len(single_allele_sols[i])
            if newmax > maxnum:
                maxnum = newmax
        logging.debug('max num of solutions: %d', maxnum)
        # branch & bound
        merged_sol = []
        numlist = [0, 1]
        merged_sol = self._merge_single_allele_solutions_lists(single_allele_sols[0], single_allele_sols[1], numlist, family.childs)
        for i, single_sol in enumerate(single_allele_sols):
            if i in (0, 1):
                continue
            numlist.append(i)
            merged_sol = self._merge_single_allele_solutions_lists(merged_sol, single_sol, numlist, family.childs)
        merged_sol = self._remove_permutations_from_states(merged_sol)
        if not merged_sol:
            logging.debug('no solution exist after branching - continue to the next family')
            self._statistics.nulled_in_merge += 1
            return None
        # DEBUG
        # logging.debug('final solutions before max len: %s', str(merged_sol))

        # need to remove permutations again cause merge non-permutations may end up in permutations of max solutions
        final_solutions = self._get_max_solutions_from_merge(merged_sol, haplotypes)
        if final_solutions is None:
            logging.debug('no solution with valid haplotypes - continue to the next family')
            self._statistics.nulled_in_max += 1
            return None
        final_solutions = self._remove_permutations_from_states(final_solutions)
        # DEBUG: print solutions
        # logging.debug(final_solutions)
        self._statistics.solved += 1
        return final_solutions

    @staticmethod
    def _remove_permutations_from_states(solutions):
        for i, sol1 in enumerate(solutions):
            if sol1 is None:
                continue
            for j, sol2 in enumerate(solutions):
                # same solutions / solution already removed
                if i == j or sol2 is None:
                    continue
                # check permutations
                if sol1.check_permutation_state(sol2):
                    # mark to remove
                    solutions[j] = None
        # remove permutations
        solutions = [x for x in solutions if x is not None]
        return solutions

    # TODO: fix for the case of unordered haplotypes/remove this
    @staticmethod
    def _fast_haplotype_check(childs):
        """
        Check if the number of different haplotypes (pairwise disjoint) in children is at most 4
        if not then there is no solution, otherwise solution maybe exist but not necessary

        Args:
            childs (list of HumanSample): children haplotypes

        """

        haplos = []
        for child in childs:
            for i in range(2):
                exists = False
                for hap in haplos:
                    if not hap.check_haplotype_disjoint(child.haplos[i]):
                        exists = True
                        break
                if not exists:
                    haplos.append(child.haplos[i])
        if len(haplos) > 4:
            return False
        return True

    @staticmethod
    def _get_not_match_alleles_list(state, child, allele_num_list):
        # intialize notmatchlist
        # notmatchlist[parent][haplo][i] = where parent's haplotype i alleles not matching to child haplotype haplo
        notmatchlist = dict()
        notmatchlist[ParentType.FATHER] = [[], []]
        notmatchlist[ParentType.MOTHER] = [[], []]
        for haplo in range(2):
            notmatchlist[ParentType.FATHER][haplo] = [[], []]
            notmatchlist[ParentType.MOTHER][haplo] = [[], []]
            for ind in range(2):
                notmatchlist[ParentType.FATHER][haplo][ind] = [[True, True]] * len(allele_num_list)
                notmatchlist[ParentType.MOTHER][haplo][ind] = [[True, True]] * len(allele_num_list)

        for ptype in [ParentType.FATHER, ParentType.MOTHER]:
            for chaplo in range(2):
                for par_hap_num in range(2):
                    major_allele_list = state.get_allele(ptype, AlleleType.MAJOR, par_hap_num, range(len(allele_num_list)))
                    minor_allele_list = state.get_allele(ptype, AlleleType.MINOR, par_hap_num, range(len(allele_num_list)))
                    if not major_allele_list or not minor_allele_list:
                        raise RuntimeError()
                    # notmatchlist[ptype][chaplo][par_hap_num] is a list contains list of pairs of True/False for example: [[True,False],[False,False],...] where the first is true if the major allele at that index covered by allele of the child and the second is for the minor
                    # this check if the child actually COVERING the parent (if the child cover the parent then the parent must match the child, the converse is not true cause the child constraint the solution)
                    notmatchlist[ptype][chaplo][par_hap_num] = child.haplos[chaplo].find_where_alleles_not_match(major_allele_list, minor_allele_list, allele_num_list)
        return notmatchlist

    @staticmethod
    def _get_match_alleles(state, child, allele_ind_list):
        """
        Check where the current solution state matches the child and return dictionary match[ptype][hapnum],
        where ptype is the parent type and hap num is the munber of the haplotype (starts from 0), if match is True then the parent allele matches the child at thtat haplotype.
        ASUUMING order of the child haplotype, child.haplos[0] = first haplotype ordered etc
        :param state: current solution state to check
        :type state: AlleleSolutionState
        :param child: child haplotype
        :type child: HaploSample
        :param allele_ind_list: allele index list in child to match
        :type allele_ind_list: list of int
        """
        # init lists of match for each haplotype
        # match[parent][i][j] is true if child's haplotype i is covered by parent haplotype j
        match = dict()
        for ptype in ParentType:
            match[ptype] = [[False, False], [False, False]]

        for ptype in ParentType:
            for ch_haplo in range(2):
                for par_haplo in range(2):
                    major = state.get_allele(ptype, AlleleType.MAJOR, par_haplo, range(len(allele_ind_list)))
                    minor = state.get_allele(ptype, AlleleType.MINOR, par_haplo, range(len(allele_ind_list)))
                    match[ptype][ch_haplo][par_haplo] = child.haplos[ch_haplo].match_alleles_list(major, minor, allele_ind_list)

        # matchparent[mother/father][i] = if mother/father alleles are matching child's allele at haplotype i
        matchparent = dict()
        matchparent[ParentType.FATHER] = [match[ParentType.FATHER][0][0] or match[ParentType.FATHER][0][1], match[ParentType.FATHER][1][0] or match[ParentType.FATHER][1][1]]
        matchparent[ParentType.MOTHER] = [match[ParentType.MOTHER][0][0] or match[ParentType.MOTHER][0][1], match[ParentType.MOTHER][1][0] or match[ParentType.MOTHER][1][1]]

        # cover nulled child haplotypes
        for hapnum in range(2):
            if all(child.haplos[hapnum].majorli[allele_ind] == '' for allele_ind in allele_ind_list):
                for ptype in ParentType:
                    matchparent[ptype][hapnum] = True

        return matchparent

    @staticmethod
    def _check_single_allele_solution_in_set(state, sol_stack):
        removeind = []
        already_in_stack = False
        for ind, sol_state in enumerate(sol_stack):
            result_state = state.compare_state(sol_state)
            # if current solution s is a subset of state
            if result_state is sol_state:
                # then s is more constrained solution of state - should replace solution s with state
                # mark the removed with None for efficient removing
                sol_stack[ind] = None
                removeind.append(ind)
            # if state is a subset of current solution s
            if result_state is state:
                # state already included in the solutions in a more constrained way (or equal), so it's already in the stack
                already_in_stack = True
                # no need to continue removing other solutions cause they already removed when result_state inserted
                break

        # efficient way to remove from list by a list of indexes to remove, by swapping
        finallen = len(sol_stack) - len(removeind)
        for i in removeind:
            if i >= finallen:
                continue
            temp = sol_stack.pop()
            while temp is None:
                temp = sol_stack.pop()
            sol_stack[i] = temp
        # now all the None need to be at the end of the sol_stack
        while len(sol_stack) != finallen:
            # temp should be None and removed
            temp = sol_stack.pop()
            if temp is not None:
                raise 'problem with None'
        if already_in_stack:
            return True
        return False

    @staticmethod
    def _check_covered(match):
        """
        Check if haplotype covered from match imformation
        :param match: describes the matches
        :type match: dictionary in the format of GetMatchSingleAllele
        """
        if (match[ParentType.MOTHER][0] and match[ParentType.FATHER][1]) or (match[ParentType.MOTHER][1] and match[ParentType.FATHER][0]):
            # current child match - continue checking the others
            return True
        return False

    @staticmethod
    # TODO: check this function
    def _check_child_covered(state, child, numlist):
        for par_order in [[0, 0], [0, 1], [1, 0], [1, 1]]:
            major1 = state.get_allele(ParentType.FATHER, AlleleType.MAJOR, par_order[0], range(len(numlist)))
            minor1 = state.get_allele(ParentType.FATHER, AlleleType.MINOR, par_order[0], range(len(numlist)))
            major2 = state.get_allele(ParentType.MOTHER, AlleleType.MAJOR, par_order[1], range(len(numlist)))
            minor2 = state.get_allele(ParentType.MOTHER, AlleleType.MINOR, par_order[1], range(len(numlist)))
            if child.match_unordered_haplos(major1, minor1, major2, minor2, numlist):
                return True
        return False

    @staticmethod
    def _get_iteration_order(fam, allele_ind):
        """
        Return the order for iteration as list.

        :param fam: family
        :type fam: Family
        """
        # iter_order save the order of the haplotypes
        iter_order = []

        # determining the iteration order for optimizations
        for parenthap in fam.father.haplos + fam.mother.haplos:
            if parenthap.minorli[allele_ind][0] != '':  # if the list not "nulled"
                # add to the start
                iter_order.insert(0, parenthap)
            else:
                # append at the end if nulled
                iter_order.append(parenthap)
        return iter_order

    @staticmethod
    def _add_next_match_state(state, match, allele_ind, child, haplo_order, method, state_stack):
        """
        Add all the next alleles in the state that cover specific child (by adding new ones or updating existing solutions) to the stack.
        :param state: current state to match
        :type state: AlleleSolutionState
        :param match: match from _check_match_single_allele
        :type match:
        :param allele_ind: the allele index that currently solved
        :type allele_ind: int
        :param child: child haplotype
        :type child: HumanSample
        :param haplo_order: the order of haplotypes returned from _get_iteration_order
        :type haplo_order:
        :param method: the method to use (update or add_allele from AlleleSolutionState)
        :type method:
        :param state_stack: the stack of the states
        :type state_stack: list of AlleleSolutionState
        """

        # if allele is added then the index should be the next order index
        if method == AlleleSolutionState.add_allele:
            ind = state.get_order_index()
        # if allele is updated then the index should be the last order index
        else:
            ind = state.get_order_index() - 1
        cur_parent_type = haplo_order[ind].haptype
        cur_hap_num = haplo_order[ind].hapnum
        # existing allele
        if method == AlleleSolutionState.update_allele:
            cur_major = state.get_allele(cur_parent_type, AlleleType.MAJOR, cur_hap_num)[0]
            cur_minor = state.get_allele(cur_parent_type, AlleleType.MINOR, cur_hap_num)[0]
        # new allele
        else:
            cur_major = haplo_order[ind].majorli[allele_ind]
            cur_minor = haplo_order[ind].minorli[allele_ind]
        for haplonum in range(2):
            # check if father major could match child major, then of course they can match
            child_major = child.haplos[haplonum].majorli[allele_ind]
            child_minor = child.haplos[haplonum].minorli[allele_ind]
            if not match[cur_parent_type][haplonum] and (cur_major in ('', child_major)):
                new_major = child_major
                new_minor = cur_minor
                if child_minor != ['']:
                    # on this case it could be only that the minor allele is free, if it's not then it should be covered by other child haplotype
                    # need to update that so the father/mother will get the intersection with the child
                    if cur_minor == ['']:
                        new_minor = child_minor[:]
                    # intersection if parent allele not free, that will give the least constrained minor allele option
                    else:
                        new_minor = list(set(cur_minor) & set(child_minor))
                    # no intersection, cant cover
                    if not new_minor:
                        continue
                newstate = copy.copy(state)
                method(newstate, new_major, new_minor, cur_parent_type, cur_hap_num)
                if method == AlleleSolutionState.add_allele:
                    newstate.increase_order_index()
                state_stack.add(newstate)

    @staticmethod
    def _solve_match_in_alleles(fam):
        """
        solves the problem of single allele and return all the solution as AlleleSolution list
        if some of the matches is empty return None

        :param fam: family to solve
        :type fam: Family
        """
        alleles_solutions = []
        alleles_num = fam.allele_count
        for allele_ind in range(alleles_num):
            sol_stack = FamilySolver._solve_single_allele_match(fam, allele_ind)
            # no solution for specific allele - therefore no solution exist
            if not sol_stack:
                return None
            alleles_solutions.append(sol_stack)
        return alleles_solutions

    @staticmethod
    def _single_allele_process_child(child, state, haplo_order, allele_ind, state_stack):
        """
        Process child covering check for single allele solution process, Returns true if child is covered by the state.
        if not covered, adding to the stack all the other possibilities for covering.
        Args:
            child (HumanSample): child haplotype
            state (AlleleSolutionState): current state of solution
            haplo_order (list of HaploSample): order of haplotypes to iterate
            allele_ind (int): the index of the allele
            state_stack (list of AlleleSolutionState): stack of states to add new possibilities

        """
        # get all the matches of current solution state with current child
        match = FamilySolver._get_match_alleles(state, child, [allele_ind])
        if FamilySolver._check_covered(match):
            # current child match - continue checking the others
            return True
        ind = state.get_order_index()
        # adding all intersection with last parent (each time only taking intersection with last parent to prevent loops)

        if ind > 0:
            FamilySolver._add_next_match_state(state, match, allele_ind, child, haplo_order, AlleleSolutionState.update_allele, state_stack)
        # adding all the next possibilities to the queue (4 is the number of the haplotypes in the solution, therefore 4 is the maximum alleles in the solution)
        if ind < 4:
            FamilySolver._add_next_match_state(state, match, allele_ind, child, haplo_order, AlleleSolutionState.add_allele, state_stack)

            cur_major = haplo_order[ind].majorli[allele_ind]
            cur_minor = haplo_order[ind].minorli[allele_ind]
            cur_haptype = haplo_order[ind].haptype
            cur_hapnum = haplo_order[ind].hapnum
            # add free allele option
            state_stack.add(copy.copy(state).add_allele(cur_major, cur_minor, cur_haptype, cur_hapnum).increase_order_index())
        return False

    @staticmethod
    def _solve_single_allele_match(fam, allele_ind):
        """
        Solves the single allele problem, returns a list of AlleleSolutionState with all possible solutions

        Args:
            fam (Family): family to solve for
            allele_ind (int): the index of the allele

        """
        sol_stack = []
        # select the iteration order
        haplo_order = FamilySolver._get_iteration_order(fam, allele_ind)
        state_stack = set()
        closed_stack = set()
        # init search stack with first parent alleles
        state_stack.add(AlleleSolutionState.create_from_haplo(haplo_order[0], allele_ind).increase_order_index())
        # sort of dfs search
        while state_stack:
            state = state_stack.pop()
            # implemented in the hash function of AlleleSolutionState
            if state in closed_stack:
                continue
            closed_stack.add(state)
            fullmatch = True
            for child in fam.childs:
                # important: don;t change the order to avoid short-circuiting
                fullmatch = FamilySolver._single_allele_process_child(child, state, haplo_order, allele_ind, state_stack) and fullmatch
            if fullmatch:
                # still need to check and add the parent's majors
                ind = state.get_order_index()
                # iterate over the remaining alleles
                for j in range(ind, len(haplo_order)):
                    state.add_allele(haplo_order[j].majorli[allele_ind], haplo_order[j].minorli[allele_ind], haplo_order[j].haptype, haplo_order[j].hapnum).increase_order_index()
                if not FamilySolver._check_single_allele_solution_in_set(state, sol_stack):
                    sol_stack.append(state)
        return sol_stack

    @staticmethod
    def _merge_single_allele_solutions_lists(sol_list1, sol_list2, numlist, childs):
        """
        Merge two allele solutions (each with different alleles solved already) into consistent solutions of the combined alleles

        Args:
            sol_list1 (list of AlleleSolutionState): list of solutions for some alleles
            sol_list2 (list of AlleleSolutionState): list of solutions for other alleles
            numlist (list of int): describing the indexes of the alleles in family haplotypes
            childs (list of HumanSample): children haplotypes

        """
        result = []
        for sol1 in sol_list1:
            for sol2 in sol_list2:
                mergedsols = FamilySolver._merge_single_allele_solutions(sol1, sol2, numlist, childs)
                for msol in mergedsols:
                    if not FamilySolver._check_single_allele_solution_in_set(msol, result):
                        result.append(msol)
        return result

    @staticmethod
    def _get_max_solutions_from_merge(merged_sol, haplotypes):
        final_max_solutions = []

        # load merged solutions into frequency haplotypes format
        for sol in merged_sol:
            final_max_solutions.append(MaxSolutionState(sol))

        # TODO: not very efficient cause haplotypes already sorted, so there is no need to iterate after solution found
        # calculates frequency and the final haplotypes with that frequency
        for hap in haplotypes:
            for msol in final_max_solutions:
                msol.process_haplo(hap)

        # now we have the maximum solutions for each merge option
        # we need to calculate the one that has maximum probability
        max_freq = 0
        # the indexes of the maximum solutions
        max_indexes = []
        for ind, msol in enumerate(final_max_solutions):
            sol_freq = msol.get_max_frequency()
            # one of the haplotypes doesn't exist or with lower frequency - remove it cause it's not a solution
            if sol_freq == 0 or sol_freq < max_freq:
                continue
            if sol_freq > max_freq:
                max_freq = sol_freq
                max_indexes = [ind]
            elif sol_freq == max_freq:
                max_indexes.append(ind)

        # no valid solution exist
        if not max_indexes:
            return None
        return [final_max_solutions[x] for x in max_indexes]

    # TODO: Split function to other functions
    @staticmethod
    def _merge_single_allele_solutions(sol1, sol2, numlist, childs):
        solutions = []
        # first insert to the stack all proposed solution from the cartesian product ignoring order
        solutions_stack = deque(sol1.add_ignoring_order(sol2))
        while solutions_stack:
            cursol = solutions_stack.popleft()
            fullmatch = True
            for child in childs:
                # check if the child already covered by the solution
                if FamilySolver._check_child_covered(cursol, child, numlist):
                    # current child match - continue checking the others
                    continue
                # one child doesn't match therefore no full match
                fullmatch = False
                # This not_match_list later used for our unordered haplotypes
                # notmatchlist[ptype][child_hap][par_hap][num][0/1] = True if ptype parent and haplotype par_hap match the child in haplo child_hap at allele number num with type 0(major) or 1(minor)
                notmatchlist = FamilySolver._get_not_match_alleles_list(cursol, child, numlist)
                # here need to copy current solution and add to the stack all 8 proposed solutions * child's ambiguity
                # for each child and order of match, there is only one way to match the parents to the child that is the least constrained.
                # iterate over all options to cover current child haplotypes - one of these should match
                for parent_order in [[0, 0], [0, 1], [1, 0], [1, 1]]:
                    # for each order need to check if parentorder[0],haploorder[0] cover the first haplotype in child and parentorder[1],haploorder[1] cover the second one
                    # parent_order[i] is the choosing for the ith parent, where i=0 is the father and i=1 is the mother.
                    # copies current state for the new solution

                    parents = [ParentType.FATHER, ParentType.MOTHER]

                    # child can have many solutions because there is no order, max of 2^ALLELES_NUM, but usually a lot (sub exponentially) less
                    # we initialize with the cursol
                    child_solutions_stack = deque([cursol])
                    # compare each allele
                    # enumerate over the alleles, choosing here the notmatch list as iterating object that doesn't really used
                    for num, ismatch in enumerate(notmatchlist[parents[0]][0][parent_order[0]]):
                        # num is the number of allele in the solution and numlist[num] is the allele number in haplosample, ismatch is a pair [majormatch,minormatch]

                        # combines the match of the major and the minor for the specific choosing of parents
                        # notmatchcombined[ptype][i] = True if the child i'th haplo (there is the first and the second each represents different haplotype) covered by the parent ptype, major+minor
                        notmatchcombined = dict()
                        for i in range(2):
                            notmatchcombined[parents[i]] = [None, None]
                            for haplo_num in range(2):
                                notmatchcombined[parents[i]][haplo_num] = notmatchlist[parents[i]][haplo_num][parent_order[i]][num][0] and notmatchlist[parents[i]][haplo_num][parent_order[i]][num][1]
                        # does the child major and minor match?
                        #logging.debug('not match for %s and allele number %s : %s', str(parent_order), str(num), str(notmatchcombined))
                        already_cover = (notmatchcombined[ParentType.FATHER][0] and notmatchcombined[ParentType.MOTHER][1]) or (notmatchcombined[ParentType.MOTHER][0] and notmatchcombined[ParentType.FATHER][1])
                        if already_cover:
                            # allele already covered cause or father cover the first and mother the second or the oposite
                            continue

                        # when testing options there are many options the could be to match the parents, but the total count cant be too high
                        # therefore we need to make a stack of all the solutions and solve them, notice that one allele solution don't change other allele solutions
                        # so what we do is to update all the solutions so far in the already processed alleles with the

                        # here we have 2 options, or father cover the first haplotype or mother cover the first and the other covered by the other parent
                        # so we have a branch factor of 2 at max from this

                        new_child_solutions_stack = deque()

                        for hap_order_num in [[0, 1], [1, 0]]:

                            # the majors and minors of the specific alleles n the parents that required for this match
                            majors = [None, None]
                            minors = [None, None]
                            # is this specific allele has no solution already?
                            curr_allele_empty = False
                            # check the cover of both parents, the 2 must have covering (this is a logical AND)
                            for i in range(2):
                                # the iteration is not depend on the other iterations, if one of solutions can't cover then no one can cover, cause every solution generated s.t it matches one match and dont change results in the other matches
                                ch_major = child.haplos[hap_order_num[i]].majorli[numlist[num]]
                                ch_minor = child.haplos[hap_order_num[i]].minorli[numlist[num]]
                                # major and minor of the parent allele that should match this child
                                major_al = cursol.get_allele(parents[i], AlleleType.MAJOR, parent_order[i], [num])[0]
                                minor_al = cursol.get_allele(parents[i], AlleleType.MINOR, parent_order[i], [num])[0]

                                # this code take the intersection in the corresponding allele
                                # major not match
                                mismatch_major = not notmatchlist[parents[i]][hap_order_num[i]][parent_order[i]][num][0]
                                # minor not match
                                mismatch_minor = not notmatchlist[parents[i]][hap_order_num[i]][parent_order[i]][num][1]
                                if mismatch_major:
                                    if major_al != '':
                                        curr_allele_empty = True
                                        break
                                    else:
                                        majors[i] = ch_major
                                else:
                                    majors[i] = major_al  # parent allele if the parent in subset of the child
                                if mismatch_minor:
                                    # if minor allele is not free and it's the only one, then the allele cannot be constrained to match
                                    if minor_al != [''] and len(minor_al) == 1:
                                        curr_allele_empty = True
                                        break
                                    # only if parent minor is [''] or a list with some alleles there could be an intersection
                                    else:
                                        # child minor can't be [''], because then there is a match in the minor
                                        # if minor_al==['']
                                        newminor = ch_minor
                                        if minor_al != ['']:
                                            # parent minor should be contained in child minor, therefore doing intersection
                                            # the sorted here is important to not over counting solutions
                                            newminor = sorted(set(ch_minor) & set(minor_al))
                                        if not newminor:
                                            curr_allele_empty = True
                                            break
                                        minors[i] = newminor
                                else:
                                    minors[i] = minor_al  # parent allele if the parent in subset of the child
                            # current allele cant match with this child hap order, skip to try the next one
                            if curr_allele_empty:
                                continue

                            # DEBUG
                            # print(majors)
                            # print(minors)
                            # print(parent_order)
                            # print(num)

                            # if we passed the loop - we found a match for the corresponding allele, so this allele solvable
                            # add the new solution to the stack - for every existing solution, also every solution usually different, complete check will occur later
                            for oldsol in child_solutions_stack:
                                sol = copy.copy(oldsol)
                                # update the alleles to match
                                for i in range(2):
                                    sol.update_allele(majors[i], minors[i], parents[i], parent_order[i], num)
                                new_child_solutions_stack.append(sol)
                        # after check all matches of child with this parent order, now new_child_solutions_stack is exactly containing the new solutions for this allele

                        # if there is some solutions update the stack and continue with them
                        child_solutions_stack = new_child_solutions_stack
                        # print(len(new_child_solutions_stack))

                        # in this allele there is no match, continue to the next parent order
                        # if one of the alleles not match - skip
                        if not child_solutions_stack:
                            break
                    # again, no solution for this parent order, skip
                    if not child_solutions_stack:
                        continue
                    # after we checked all alleles we count the new solution for this parent
                    # DEBUG
                    #print('child solutions= ', len(child_solutions_stack))
                    # here we have solutions for this parent order, this solution cover at least 1 child, add that to the solutions_stack
                    for newsol in child_solutions_stack:
                        if not FamilySolver._check_single_allele_solution_in_set(newsol, solutions_stack):
                            solutions_stack.append(newsol)
                # DEBUG
                # print('solutions len = ', len(solutions_stack))
                # print('solutions = ', solutions_stack)
            if fullmatch:
                solutions.append(cursol)
        return solutions
