'''
Contains serialized file ambiguity reader class
'''

import os.path
import pickle
import array
from reader.ambiguity_reader import AmbiguityReader


class SerializedFileAmbiguityReader(AmbiguityReader):
    """
    Ambiguity reader implementation for serialized files
    """

    SERIALIZATION_EXT = 'ser'

    def __init__(self, full_filename):
        """
        Creates reader using the file
        :param full_filename: path to the file
        :type full_filename: string
        """

        filename, _ = os.path.splitext(full_filename)
        self._serialized = False
        self._serial_filename = filename + '.' + self.SERIALIZATION_EXT
        if os.path.isfile(self._serial_filename):
            self._serialized = True
            data = open(self._serial_filename, 'rb')
        else:
            self._serialized = False
            data = open(full_filename, 'r')
        super().__init__(data)

    def _get_next_entry(self):
        """
        if serialized then automatically load the table and return None
        """
        if self._serialized:
            self.translation = pickle.load(self._data)
            return None
        return next(self._data, None)

    def _get_ambiguity_key(self, entry_data):
        if self._serialized:
            raise RuntimeError()
        return entry_data.split(' ')[0]

    def _get_ambiguity_alleles(self, entry_data):
        if self._serialized:
            raise RuntimeError()
        # the second part of the contain the alleles
        alleles_str = entry_data.split(' ')[1].rstrip('\n')
        # the different options in the alleles
        options = alleles_str.split('/')
        # if the major allele defines the ambiguity
        if ':' in alleles_str:
            dict_al_maj = dict()
            # now options is a list of [major,minor] for the different alleles
            options = [[x.split(':')[0], x.split(':')[1]] for x in options]
            for option in options:
                # major allele option not added yet
                if int(option[0]) not in dict_al_maj:
                    dict_al_maj[int(option[0])] = array.array('H')
                # remove last character from minor (17q -> 17)
                if not option[1].isdigit():
                    option[1] = option[1][:-1]

                dict_al_maj[int(option[0])].append(int(option[1]))
            return dict_al_maj
        # ambiguity independent in major allele
        for i, option in enumerate(options):
            # remove last character from minor (17q -> 17)
            if not options[i].isdigit():
                options[i] = options[i][:-1]
        # array of shorts saves a lot of memory
        return array.array('H', [int(x) for x in options])

    def _free_resources(self):
        """
        Also used to save serialization data file
        """
        if not self._serialized:
            serialfile = open(self._serial_filename, 'wb')
            pickle.dump(self.translation, serialfile)
            serialfile.close()
        self._data.close()
