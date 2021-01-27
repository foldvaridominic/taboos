import logging
import random
from collections import Counter, defaultdict
from functools import reduce
from itertools import combinations


logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
ch.setLevel(logging.INFO)
logger.addHandler(ch)


from constants import NUCLEOTIDE_CHARACTERS, ALPHABET_LENGTH
from tabooset import TabooSet as TS
from utils import direct_product, to_string


letters = set(NUCLEOTIDE_CHARACTERS)
max_letters = ALPHABET_LENGTH
allowed = []

def rec_func(length):
    global allowed
    for r in range(1,max_letters + 1):
        allowed.append(set(random.sample(letters, r)))
        if length > 1:
            yield from rec_func(length-1)
        else:
            yield allowed
        allowed = allowed[:-1]


def enum_graphs(length):
    all_blocks = set()
    all_structures = defaultdict(list)
    for allowed in rec_func(length):
        not_allowed = [letters-pos for pos in allowed]
        taboos = []
        for idx in range(length):
            taboos += [to_string(s) for s in direct_product(allowed[:idx] + [not_allowed[idx]] + allowed[idx+1:])]
        structure = tuple([len(a) for a in allowed])
        logger.info("Taboo blocks: %s", structure)
        blocks = TS.gen_hamming_graph(taboos, length)
        if len(blocks) > 1:
            block = frozenset(Counter(blocks).items())
            all_blocks.add(block)
            all_structures[block].append(structure)
    logger.info("All different block structures: %s", all_blocks)
    logger.info("Number of different block structures: %s", len(all_blocks))
    #logger.info(all_structures)


def sigma_minus(number, alphabet_length=4):
    return alphabet_length - number


def flip_at(iterator, indices):
    ret = [i for i in iterator]
    for i in indices:
        ret[i] = sigma_minus(ret[i])
    return ret


def enum_quotient_graph(dimensions, alphabet_length=4):
    length = len(dimensions)
    node_levels = []
    for flip in range(length+1):
        node_levels.append(sum(reduce(lambda x,y: x*y,  flip_at(dimensions, indices)) for indices in combinations(range(length),flip)))
    return node_levels


def increase_dimension(characters, increase=4):
    not_allowed = [letters-pos for pos in characters]
    taboos = []
    for idx, _ in enumerate(characters):
        taboos += [to_string(s) for s in direct_product(characters[:idx] + [not_allowed[idx]] + characters[idx+1:])]
    structure = tuple([len(c) for c in characters])
    logger.info("Taboo blocks: %s", structure)
    for jdx in range(increase):
        blocks = TS.gen_hamming_graph(taboos, idx+1+jdx)
    return
