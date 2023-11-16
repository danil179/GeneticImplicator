# pylint: disable=C0301

'''
Contains table data file reader class.
'''

import logging
from openpyxl import load_workbook
from reader.data_reader import DataReader
from genetic.human_sample import HumanSample
from genetic.haplo_sample import HaploSample
from genetic.family import Family
from genetic.parent_type import ParentType
from genetic.allele_type import AlleleType


class TableFileDataReader(DataReader):
    """Read the genetics data from table file"""

    def __init__(self, file_path, translation_table, keep_nulled_childs=False, remove_redundant_childs=False):
        """
        Creates a file reader that reads from table

        :param file_path: path to the file used as 'data' in DataReader
        :type file_path:
        :param translation_table: dictionary that take ambiguity string and returns a list of integers (the alleles)
        :type translation_table: TranslationTable
        :param nulled_childs:
        :type nulled_childs:
        :param redundant_childs:
        :type redundant_childs:
        """
        super().__init__(file_path, keep_nulled_childs, remove_redundant_childs)
        self._init_table_data()
        self._transtable = translation_table

    def _init_table_data(self):
        """Initiate table data and variables"""
        self._wb = load_workbook(filename=self.data, read_only=True)
        self._ws = self._wb['Sheet1']
        self._fam_code_name = 'FAMCODE'
        self._type_name = 'BIRTHSEQ'
        self._major_allele_names = [['A1', 'B1', 'C1', 'DR1', 'DQ1'], ['A2', 'B2', 'C2', 'DR2', 'DQ2']]
        self._minor_allele_names = [['HLA-A*1', 'B_STAR1', 'C_STAR1', 'DRB11', 'DQB11'], ['HLA-A*2', 'B_STAR2', 'C_STAR2', 'DRB12', 'DQB12']]
        self.__col_names = ['FAMCODE', 'SEQNO', 'BIRTHSEQ', 'HLA-A*1', 'HLA-A*2', 'A1', 'A2', 'B_STAR1', 'B_STAR2', 'B1', 'B2', 'C_STAR1', 'C_STAR2', 'C1', 'C2', 'DRB11', 'DRB12', 'DR1', 'DR2', 'DQB11', 'DQB12', 'DQ1', 'DQ2']
        self._cur_row = 0
        self.fam_count = 0
        self._rows = list(self._ws.rows)
        self._max_row = len(self._rows)
        header = self._rows[0]
        self._cur_row += 1
        self.__col_indices = {cell.value: n for n, cell in enumerate(header) if cell.value in self.__col_names}

    @staticmethod
    def __read_major(major_str):
        """
        This function reads the string of the major part and return it as int or '' (for blank)

        :param major_str: major string (e.g. 17q)
        :type major_str: string
        """
        if major_str is None or major_str == '':
            return ''
        # remove last character like p/n/q
        if not major_str.isdigit():
            if not major_str[:-1].isdigit() or len(major_str) < 2:
                logging.info('unknown major symbol: %s - considered blank', major_str)
                return ''
            major_str = major_str[:-1]
        return int(major_str)

    @staticmethod
    def __read_major_part_of_minor(major_part_str):
        """
        Converts major part of minor and returns it as int or '' (for blank)

        Args:
            major_part_str (string): major string (e.g. 17q)

        """
        if major_part_str is None or major_part_str == '':
            return ''
        # remove last character like p/n/q if the minor part
        if not major_part_str.isdigit():
            major_part_str = major_part_str[:-1]
        return int(major_part_str)

    def __read_minor_part_of_minor(self, minor_part_str, major_allele):
        """
        Converts minor part of minor and return it as list of int or None (for unknown ambiguity)

        :param minor_part_str: minor string (e.g. 17,AB)
        :param major_allele: major allele
        """

        if minor_part_str[0].isupper():  # if minor part starts with uppercase letter then there is ambiguity (e.g. AB)
            return self._transtable.resolve_minor(minor_part_str, major_allele)
        # minor part of the allele is known (e.g. 14p)
        if not minor_part_str.isdigit():
            minor_part_str = minor_part_str[:-1]
        return [int(minor_part_str)]

    @staticmethod
    def __get_major_allele_consistent(majorl1, majorl2):
        """
        Checks if alleles consistent with the rules of inference from this data type
        Args:
            majorl1 (list of (int/'')): one list of majors
            allele2 (list of (int/'')): one list of majors
        """

        permutation_order = None
        # check permutation equality
        for order in [[0, 1], [1, 0]]:
            skip = False
            for order1, order2 in enumerate(order):
                if TableFileDataReader.__intersect_allele(majorl1[order1], majorl2[order2]) is None:
                    # no intersection - check next permutation
                    skip = True
                    break
            if skip:
                continue
            # found permutation - the other permutation could give another options, but when we duplicate alleles this is not the case
            permutation_order = order
            break

        if permutation_order is None:
            # inconsistent
            return None

        major = [TableFileDataReader.__intersect_allele(majorl1[i], majorl2[permutation_order[i]]) for i in range(2)]

        return major

    @staticmethod
    def __intersect_allele(al1, al2):
        """
        Intersection between alleles, returns the intersection if not empty, else returns None.
        Args:
            major1 (int/''): major allele 1
            major2 (int/''): major allele 2
        """

        if al1 == '':
            return al2
        if al2 == '':
            return al1

        if al1 != al2:
            return None

        return al1

    def _read_row(self, row):
        '''
        Read row from the table into 2 haplotypes dictionary, if there is some error return None.

        :param row: the row to convert to haplotypes
        :type row: list of Cell, returned from rows[index] with openpyl
        '''
        # haplo contains the different haplotypes
        haplo = dict()
        haplo[AlleleType.MAJOR] = [[], []]
        haplo[AlleleType.MINOR] = [[], []]
        # TODO: configurable column names
        # column names in the table
        for allele_ind in range(len(self._major_allele_names[0])):
            # read major allele column, major[i] is the allele from the ith haplo
            major_allele = [self.__read_major(row[self.__col_indices[self._major_allele_names[i][allele_ind]]].value) for i in range(2)]
            # read minor allele column
            minor_details = [row[self.__col_indices[self._minor_allele_names[i][allele_ind]]].value for i in range(2)]
            # replace None cells with ''
            minor_details = ['' if minor_details[i] is None else minor_details[i] for i in range(2)]
            # spliting the major and minor parts, e.g. minor_parts = [[12,1],[13,2]]
            minor_parts = [minor_details[i].split(':') for i in range(2)]
            # read major's part in the minor column
            major_allele_minor_version = [self.__read_major_part_of_minor(minor_parts[i][0]) for i in range(2)]
            # checks if major allele consistent from the other parts
            major = self.__get_major_allele_consistent(major_allele, major_allele_minor_version)
            if major is None:
                # inconsistent
                return None
            # now major is a list of the majors, the order could be different from the minors

            # init minor as size 2 list for each haplotype
            minor = [None] * 2
            for i in range(2):
                # no minor part
                if len(minor_parts[i]) == 1:
                    minor[i] = ['']
                # minor part exists
                else:
                    # major_allele_minor_version used here for the correct order
                    minor[i] = self.__read_minor_part_of_minor(minor_parts[i][1], major_allele_minor_version[i])
                # error in minor part format
                if minor[i] is None:
                    return None
            # now minor is a list of the minors, the order could be different from the majors

            # finding the order, if two minors = [''] the order is not important
            final_order = [0, 1]
            for i in range(2):
                if minor[i] != ['']:
                    # major_allele_minor_version must be some allele, compraing order
                    if major_allele_minor_version[i] == major[i]:
                        final_order = [0, 1]
                    # the order is swapped
                    else:
                        final_order = [1, 0]

            # duplicating alleles
            for i in range(2):
                if major[i] == '':
                    j = 1 if i == 0 else 0
                    # two major alleles empty
                    if major[j] == '':
                        break
                    # only this allele is empty - duplicating the other to this allele
                    major[i] = major[j]
                    # important that we duplicate inorder
                    minor[final_order[i]] = minor[final_order[j]]
                elif minor[final_order[i]] == ['']:
                    j = 1 if i == 0 else 0
                    if major[j] != major[i] or minor[final_order[j]] == ['']:
                        continue
                        # two major alleles equal and the other minor not empty - assumes that both minor alleles the same
                    minor[final_order[i]] = minor[final_order[j]]

            # no order preserving in haplotypes
            for i in range(2):
                haplo[AlleleType.MAJOR][i].append(major[i])
                haplo[AlleleType.MINOR][i].append(minor[final_order[i]])
        # haplo[alleletype (major/minor)][haplotype index] = haplotype major/minor alleles with corresponding index
        return haplo

    def read_next_family(self):
        """Read the next family using some default format and return Family"""
        # end of table reached
        if self._cur_row >= self._max_row:
            return None
        # read family code
        row = self._rows[self._cur_row]
        fam_code = row[self.__col_indices['FAMCODE']].value
        # create family container
        fam = Family(fam_code)

        # do-while loop - check end of table not reached
        while self._cur_row < self._max_row:
            row = self._rows[self._cur_row]
            new_fam_code = row[self.__col_indices['FAMCODE']].value
            # new family starts - breaks
            if new_fam_code != fam.code:
                break
            self._cur_row += 1

            haplo = self._read_row(row)
            if not haplo:
                logging.info('inconsistent line %d', self._cur_row)
                continue

            # t is the type of the haplotype ('F','M','C')
            ptype = row[self.__col_indices['BIRTHSEQ']].value
            seq_num = row[self.__col_indices['SEQNO']].value
            # add the human to the family
            if ptype not in (ParentType.FATHER.value, ParentType.MOTHER.value):
                ptype = None
            else:
                ptype = ParentType(ptype)
            hsample = HumanSample(HaploSample(haplo[AlleleType.MAJOR][0], haplo[AlleleType.MINOR][0], ptype, 0), HaploSample(haplo[AlleleType.MAJOR][1], haplo[AlleleType.MINOR][1], ptype, 1))
            if ptype is None:
                child = hsample
                if not child.is_empty() or self.read_nulled_childs:
                    fam.add_child(child, seq_num, self.remove_redundent)
            if ptype == ParentType.FATHER:
                fam.add_father(hsample)
            if ptype == ParentType.MOTHER:
                fam.add_mother(hsample)
        # Add father and mother if none
        if fam.father is None:
            fam.add_father(HumanSample(HaploSample(['', '', '', '', ''], [[''], [''], [''], [''], ['']], ParentType.FATHER, 0), HaploSample(['', '', '', '', ''], [[''], [''], [''], [''], ['']], ParentType.FATHER, 1)))
        if fam.mother is None:
            fam.add_mother(HumanSample(HaploSample(['', '', '', '', ''], [[''], [''], [''], [''], ['']], ParentType.MOTHER, 0), HaploSample(['', '', '', '', ''], [[''], [''], [''], [''], ['']], ParentType.MOTHER, 1)))
        self.fam_count += 1
        return fam
