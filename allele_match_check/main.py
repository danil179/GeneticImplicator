'''
main module
'''

import sys
import logging
import time
import datetime
import constants
from genetic.haplo_freq_sample import HaploFreqSample
from reader.table_file_data_reader import TableFileDataReader
from reader.serialized_file_ambiguity_reader import SerializedFileAmbiguityReader
from reader.haplo_allele_reader import HaploAlleleReader
from writer.haplo_allele_writer import HaploAlleleWriter
from writer.solution_writer import SolutionWriter
from solver.family_solver import FamilySolver

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
HAPLOTYPES_FILENAME = 'haplotypes.txt'
FREQ_FILENAME = 'freq.txt'
AMBIGUITY_FILENAME = 'uncertain.txt'
DATA_FILENAME = 'DD.xlsx'
OUTPUT_FILE = 'solution.csv'

START_TIME = time.time()
FILE_HAPLO = open(HAPLOTYPES_FILENAME, 'r')
FILE_FREQ = open(FREQ_FILENAME, 'r')
HAPLO_LINES = FILE_HAPLO.read().splitlines()
FREQ_LINES = FILE_FREQ.read().splitlines()
HAPLO_ALLELE_NAME_ORDER = ['A', 'B', 'C', 'DRB1', 'DQB1']
HAPLO_ALLELE_READER = HaploAlleleReader(HAPLO_ALLELE_NAME_ORDER, '*', '~')
HAPLO_ALLELE_WRITER = HaploAlleleWriter(HAPLO_ALLELE_NAME_ORDER, '*', '~')
SOLUTION_WRITER = SolutionWriter(OUTPUT_FILE, HAPLO_ALLELE_WRITER)
HAPLOTYPES = HAPLO_ALLELE_READER.read_haplotypes(HAPLOTYPES_FILENAME)
HAPLOTYPES_PROBABILITY = []
for i, haplo in enumerate(HAPLOTYPES):
    HAPLOTYPES_PROBABILITY.append(HaploFreqSample(HAPLOTYPES[i], int(FREQ_LINES[i])))
    constants.TOTAL_FREQ += int(FREQ_LINES[i])
logging.debug('sum of all frequencies: %d', constants.TOTAL_FREQ)
# sort in decresing probability order
HAPLOTYPES_PROBABILITY = sorted(HAPLOTYPES_PROBABILITY, key=lambda l: -l.frequency)
# read ambiguities into the translation table
AMBIGUITY_READER = SerializedFileAmbiguityReader(AMBIGUITY_FILENAME)
AMBIGUITY_READER.read_ambiguities()
# get families
# TODO: add configurable table field names
DATA_READER = TableFileDataReader(DATA_FILENAME, AMBIGUITY_READER.translation, False, True)
FAMILIES = DATA_READER.read_families()
SOLVER = FamilySolver()
for fam in FAMILIES:
    # DEBUG
    # if int(fam.code) < 131037:
    #    continue
    famstime = time.time()
    logging.debug('fam code: %s', fam.code)

    solutions = SOLVER.solve(fam, HAPLOTYPES_PROBABILITY)
    if solutions is not None:
        SOLUTION_WRITER.write_solution(solutions, fam)

    logging.debug('time: %s', str(datetime.timedelta(seconds=time.time() - famstime)))
    logging.debug('finish family')
logging.debug('running time: %s', str(datetime.timedelta(seconds=time.time() - START_TIME)))
# REVIEW: consider writing stats of the solver to file
