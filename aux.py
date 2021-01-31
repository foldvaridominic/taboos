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


def enum_quotient_graph_length(dimensions, alphabet_length=4):
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


def enum_quotient_graph(characters, skip=1):
    not_allowed = [letters-pos for pos in characters]
    length = len(characters)
    node_groups = []
    index_groups = []
    jdx = 0
    for idx in range(length+1):
        if idx == skip:
            continue
        index_groups.append([])
        for indices in combinations(range(length),idx):
            jdx += 1
            index_groups[idx].append(jdx)
            current_characters = [c for c in characters]
            for i in indices:
                current_characters[i] = not_allowed[i]
            node_groups.append((jdx,current_characters))
    index_groups = [(1,), tuple(range(2,jdx+1))]
    return node_groups, index_groups


def extend_left_and_right(node_groups):
    node_groups_right = []
    node_groups_left = []
    for group in node_groups:
        node_groups_left.append((group[0], [letters] + group[1]))
        node_groups_right.append((group[0], group[1] + [letters]))
    return node_groups_left, node_groups_right


def overlaps(node_groups, index_groups):
    left, right = extend_left_and_right(node_groups)
    new_node_groups =  []
    for p1,p2 in direct_product([left, right]):
        ints = [s1.intersection(s2) for s1, s2 in zip(p1[1],p2[1])]
        if all(ints):
            new_node_groups.append(((p1[0], p2[0]), ints))
    indices = [n[0] for n in new_node_groups]
    new_index_groups = [] #TODO
    return new_node_groups, new_index_groups


def inspect_dimension_increment_in_quotient_graph(characters, inc=1):
    start_node_groups, start_index_groups = enum_quotient_graph(characters)
    node_groups, index_groups = overlaps(start_node_groups, start_index_groups)
    for s in start_node_groups:
        logger.info(s)
    for n in node_groups:
        logger.info(n)
