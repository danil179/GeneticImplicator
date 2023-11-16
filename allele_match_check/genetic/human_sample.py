# pylint: disable=C0301

'''
Contains human sample class
'''


class HumanSample:
    '''
    Describing genetic sample of human (including 2 haplotypes).
    '''

    def __init__(self, _haplo1, _haplo2):
        self.haplos = []
        self.haplos.append(_haplo1)
        self.haplos.append(_haplo2)

    @property
    def allele_count(self):
        """
        Allele count property, returns the number of alleles in the haplotypes (haplotype's allele count should be equal).
        """
        return self.haplos[0].allele_count

    def is_empty(self):
        """
        Returns true if both haplotypes empty (means that all alleles are empty string)
        """
        return self.haplos[0].is_empty() and self.haplos[1].is_empty()

    def match_allele_at(self, major, minor, i):
        """
        Return true if any haplotype match to the specific allele

        Args:
            major (int): major part of allele
            minor (int): minor part of allele
            i (type): allele index according to the haplotype
        """
        return self.haplos[0].match_allele_at(major, minor, i) or self.haplos[1].match_allele_at(major, minor, i)

    def match_unordered_haplos(self, majorl1, minorl1, majorl2, minorl2, numlist):
        """
        Return true if each major & minor in majorl & minorl are covered by some (unordered) haplotype allele (like inclusion in sets), hence the part of the haplotype match.

        Args:
            majorl (list of int): list of the allele majors
            minorl (list of list of int): list of the allele minors
            numlist (list of int): indexes to comapre with haplotype

        """
        # empty list means that no allele assigned therefore no match occur
        if majorl1 == [] or minorl1 == [] or majorl2 == [] or minorl2 == []:
            return False
        for i, num in enumerate(numlist):
            # check that there is a match with one and the other match the second
            if not ((self.haplos[0].match_allele_at(majorl1[i], minorl1[i], num) and self.haplos[1].match_allele_at(majorl2[i], minorl2[i], num))
                    or (self.haplos[1].match_allele_at(majorl1[i], minorl1[i], num) and self.haplos[0].match_allele_at(majorl2[i], minorl2[i], num))):
                return False
        return True

    def match_unordered_haplos_samp(self, haplo_samp1, haplo_samp2):
        """
        Same as match_unordered_haplos, only different parameters type

        Args:
            haplo_samp (HaploSample): haplotype sample to compare with
        """
        majorl1 = haplo_samp1.majorli
        minorl1 = haplo_samp1.minorli
        majorl2 = haplo_samp2.majorli
        minorl2 = haplo_samp2.minorli
        if (len(majorl1) != len(majorl2)) or (len(minorl1) != len(minorl2)):
            raise Exception("match_unordered_haplos - wrong parameters len")
        numlist = range(len(majorl1))

        return self.match_unordered_haplos(majorl1, minorl1, majorl2, minorl2, numlist)

    def check_sample_is_redundant(self, sample):
        '''
        returns the most constrained child if one of the children is covered by another
        '''
        for order in [[0, 1], [1, 0]]:
            # if sample covering current haplo
            if sample.haplos[order[0]].cover(self.haplos[0]) and sample.haplos[order[1]].cover(self.haplos[1]):
                return self
            # if haplo covering sample
            if self.haplos[0].cover(sample.haplos[order[0]]) and self.haplos[1].cover(sample.haplos[order[1]]):
                return sample
        # no haplosample covering another
        return None
