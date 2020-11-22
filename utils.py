import itertools
import re


def flatten(iterable):
    return itertools.chain.from_iterable(iterable)


def direct_product(iterable):
    return itertools.product(*iterable)


def get_self_product(iterable, reps):
    return itertools.product(iterable, repeat=reps)


def combinations(iterable, reps):
    return itertools.combinations(iterable, reps)


def combinations_with_replacement(iterable, reps):
    return itertools.combinations_with_replacement(iterable, reps)


def hamming_distance_1_for_tuples(pair):
    return sum(1 for t1, t2 in zip(*pair) if not set(t1).intersection(set(t2))) <= 1


def hamming_distance_1_for_strings(pair):
    return sum(1 for t1, t2 in zip(*pair) if t1 != t2) <= 1


def to_string(iterable, concat=''):
    return concat.join(iterable)
 
 
def taboo_free_for_strings(word, taboos):
    return all(not re.search(t, word) for t in taboos)


def taboo_free_for_tuples(word, taboos):
    # works as a prefilter
    return any(all(not re.search(t, to_string(p)) for t in taboos) for p in direct_product(word))


def match_from_start(pattern, word):
    return re.match(pattern, word)
