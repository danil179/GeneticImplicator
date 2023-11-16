# pylint: disable=C0301

'''
Contains haplotype sample class
'''


class HaploSample:
    '''
    Haplotype sample (describing the haplotype as list of alleles, each allele contains 2 parts: major, minor)
    '''
    HAPLO_MATCH_SIZE_ERROR = 'Cannot match haplotypes of different size'

    def __init__(self, major, haplosamp, ht, hnum):
        '''
        haplosamp is a list of lists, each list describes optional matches with the 2nd part of the corresponding allele
        major is a list of shorts for each allele
        '''
        self.majorli = major
        self.minorli = haplosamp
        self.haptype = ht
        self.hapnum = hnum

    def set_haplo_type(self, hap_type):
        """
        Set the haplotype type (if it's from father, mother or child)

        Args:
            hap_type (ParentType): haplotype type

        """
        self.haptype = hap_type

    @classmethod
    def create_from_string(cls, haplostr, reader):
        """
        Create haplotype from string using a HaploReader.

        Args:
            haplostr (string): haplotype string representation
            reader (HaploReader): convertor from string to haplotype

        """
        return reader.read_haplo_from_string(haplostr)

    def cover(self, haplo):
        """
        Return True if haplo is subset of self haplotype, else return false.

        Args:
            haplo (HaploSample): haplotype for checking covering

        """
        if len(haplo.minorli) != len(self.minorli):
            raise RuntimeError(self.HAPLO_MATCH_SIZE_ERROR)
        for i in range(haplo.allele_count):
            if self.majorli[i] != '' and self.majorli[i] != haplo.majorli[i]:
                return False
            if '' not in self.minorli[i] and not all(al in self.minorli[i] for al in haplo.minorli[i]):
                return False
        return True

    def check_haplotype_disjoint(self, haplo):
        """
        Returns true if no haplotype exist in the intersection of the haplotypes.

        Args:
            haplo (HaploSample): haplotype to check

        """
        if len(haplo.minorli) != len(self.minorli):
            raise RuntimeError(self.HAPLO_MATCH_SIZE_ERROR)
        for i in range(self.allele_count):
            if haplo.majorli[i] != '' and self.majorli[i] != '' and self.majorli[i] != haplo.majorli[i]:
                return True
            # each list contains some options and there isn't intersection
            if '' not in self.minorli[i] and '' not in haplo.minorli[i] and not list(set(self.minorli[i]) & set(haplo.minorli[i])):
                return True
        return False

    def find_where_alleles_not_match(self, majorlist, minorlist, numlist):
        """
        Returns list of pairs (as list with length of 2), each index in the list corresponds to allele in the haplotypes,
        where the first index in the pair is True if major part matches, and the second one true if minors match.

        Args:
            majorlist (list of int): major part list
            minorlist (list of list of int): minor part list
            numlist (list of int): length as majorlist, each index corresponding to allele index in the haplotype (e.g numlist[1]=2 meaning majorlist[1] will be compared with index 2 at the haplotype)

        """
        matchtable = []
        for i, num in enumerate(numlist):
            majormatch = False
            if self.majorli[num] == '' or self.majorli[num] == majorlist[i]:
                majormatch = True
            match = [majormatch, self.match_minor_allele_at(minorlist[i], num)]
            matchtable.append(match)
        return matchtable

    def match_minor_allele_at(self, minorlist, i):
        """
        Returns true if minorlist is subset of the ambiguities in index i of haplotype minor, that checks if all the options match on that minor.

        Args:
            minorlist (list of int): ambiguities ti match
            i (int): allele index according to haplotype

        """
        if '' not in self.minorli[i] and not set(minorlist).issubset(set(self.minorli[i])):
            return False
        return True

    def match_allele_at(self, major, minor, i):
        """
        Return true if the major & minor are covered by haplotype allele (like inclusion in sets), hence the allele match.

        Args:
            major (int): major part of allele
            minor (list of int): minor part of allele
            i (int): allele index to campare according to the haplotype

        """
        if self.majorli[i] == '' or self.majorli[i] == major:
            return self.match_minor_allele_at(minor, i)
        return False

    def match_alleles_list(self, majorl, minorl, numlist):
        """
        Return true if each major & minor in majorl & minorl are covered by haplotype allele (like inclusion in sets), hence the part of the haplotype match.

        Args:
            majorl (list of int): list of the allele majors
            minorl (list of list of int): list of the allele minors
            numlist (list of int): indexes to comapre with haplotype

        """
        # empty list means that no allele assigned therefore no match occur
        if majorl == [] or minorl == []:
            return False
        for i, num in enumerate(numlist):
            if not self.match_allele_at(majorl[i], minorl[i], num):
                return False
        return True

    @property
    def allele_count(self):
        """
        Allele count in the haplotype (major should be equal to minor)
        """
        return len(self.majorli)

    def __repr__(self):
        return str(self)

    def __str__(self):
        return str([self.majorli, self.minorli, self.haptype, self.hapnum])

    def is_empty(self):
        """
        Return true if haplotype empty (all the alleles are '')
        """
        for major in self.majorli:
            if major != '':
                return False
        for options in self.minorli:
            if options and '' not in options:
                return False
        return True
