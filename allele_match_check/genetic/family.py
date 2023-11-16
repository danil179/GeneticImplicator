'''
Contains family class.
'''


class Family:
    '''
    Describes family genetic structure (including all childs and parents).
    '''

    ALLELE_COUNT_ERROR = 'Allele count doesnt match the existing haplotypes'

    def __init__(self, code, _father=None, _mother=None, _childs=None):
        self.code = code
        self.father = _father
        self.mother = _mother
        self.childs = _childs
        self.full_ordered_childs = _childs
        if _childs is None:
            self.child_indexes = []
        else:
            self.child_indexes = [i for i, _ in enumerate(_childs)]
        self._allele_count = None
        if self.childs is None:
            self.childs = []
        if self.full_ordered_childs is None:
            self.full_ordered_childs = []

    def add_child(self, nchild, seq_num, rem_redundant):
        """
        Add child to the family.

        :param nchild: the new child to add
        :type nchild: HumanSample
        :param rem_redundant: check if the child redundant and remove redundancy
        :type rem_redundant: boolean
        """
        self.full_ordered_childs.append(nchild)
        self.child_indexes.append(seq_num)
        if rem_redundant:
            for i, child in enumerate(self.childs):
                rchild = child.check_sample_is_redundant(nchild)
                if rchild is not None:
                    self.childs[i] = rchild
                    return
        self.childs.append(nchild)
        self._set_allele_count(nchild)

    # TODO: replace with add_parent
    def add_father(self, nfather):
        """
        Adds father to the family if not exist yet.

        Args:
            nfather (HumanSample): father human sample

        """
        if self.father is None:
            self.father = nfather
            self._set_allele_count(self.father)

    def add_mother(self, nmother):
        """
        Adds mother to the family if not exist yet.

        Args:
            nmother (HumanSample): mother human sample

        """
        if self.mother is None:
            self.mother = nmother
            self._set_allele_count(self.mother)

    def _set_allele_count(self, haplo):
        """
        Enforcing all haplotypes to be with the same length.

        Args:
            haplo (HaploSample): haplotype to check allele length

        """
        if self.allele_count is not None and self.allele_count != haplo.allele_count:
            raise RuntimeError(self.ALLELE_COUNT_ERROR)
        self._allele_count = haplo.allele_count

    @property
    def allele_count(self):
        """
        Count of the alleles in the haplotypes (should be equal for all the haplotypes).
        """
        return self._allele_count
