import sys
import argparse
import logging
import time
import random
import itertools
import re

from constants import *


class TabooSet:
    def __init__(self, taboos, idx):
        self.idx = idx
        logger.debug("%s taboo init data: %s", idx, taboos)
        self.taboos = set(self.transform(taboos))
        #logger.debug("transformed taboos: %s", self.taboos)
        self.complement()
        self.minimize()
        self.taboo_strings = {''.join(t) for t in self.taboos}
        self.taboo_inits = {t[0] for t in self.taboos if t}
        self.taboo_init_counts = len(self.taboo_inits)

    def transform(self, taboos):
        for t in taboos:
            yield from self.generate_variants(t)

    def complement(self):
        complements = set()
        for t in self.taboos:
            complements.add(self.generate_complement(t))
        #logger.debug("complemented taboos: %s", complements)
        self.taboos = self.taboos.union(complements)

    def minimize(self):
        #TODO
        pass

    def generate_variants(self, taboo_string):
        variants = [
            (l,) if l in NUCLEOTIDE_CHARACTERS else SPECIAL_CHARACTERS.get(l)
            for l in taboo_string]
        # ignore all other characters for now
        variants_2 = [v for v in variants if v is not None]
        diff = len(variants) - len(variants_2)
        if diff:
            logger.debug("Unknown characters in %s", taboo_string)
        taboo_strings = itertools.product(*variants_2)
        yield from taboo_strings

    def generate_complement(self, taboo_string):
        return tuple(NUCLEOTIDE_COMPLEMENTS[l] for l in taboo_string[::-1])

    @property
    def connected_1(self):
        return self.taboo_init_counts < ALPHABET_LENGTH

    @property
    def connected_2(self):
        sorted_taboos = set(sorted([t[1:] for t in self.taboos], key=len, reverse=True))
        pairs = list(itertools.combinations(sorted_taboos, 2))
        #logger.debug("Taboo pairs: %s", pairs)
        pairs_filtered = [p for p in pairs if sum(1 for t1,t2 in zip(*p) if t1 != t2) <= 1]
        logger.debug("Taboo pairs filtered: %s", pairs_filtered)
        for pair in pairs_filtered:
            for idx, letter in enumerate(NUCLEOTIDE_CHARACTERS, 1):
                first = letter + ''.join(pair[0])
                second = letter + ''.join(pair[1])
                if all(t not in first and t not in second for t in self.taboo_strings):
                    break
                else:
                    logger.debug("%s or %s in taboo set", first, second)
                if idx == ALPHABET_LENGTH:
                    return False
        return True

    @property
    def connected_3(self):
        #TODO
        return False

    @property
    def connected(self):
        if self.connected_1:
            logger.debug("Connected 1: %s", self.taboo_inits)
            return True
        if self.connected_2:
            logger.info("Connected 2")
            return True
        if self.connected_3:
            return True
        return False


def run(organism_regex, enzyme_regex, rsite_regex, max_wildcard):

    idx=0
    start_time = time.time()
    re_organism = re.compile(organism_regex)
    re_enzyme = re.compile(enzyme_regex)
    re_rsite = re.compile(rsite_regex)
    organism = None
    data = []
    collect = False
    not_connected_count = 0


    for line in sys.stdin:
        name_match = re_organism.search(line)
        enzyme_match = re_enzyme.search(line)
        #XXX consider all entries as type 2 for now
        # current regex cannot be applied
        enzyme_match = True
        rsite_match = re_rsite.search(line)

        if enzyme_match:
            collect = True
        if name_match and organism is None:
            organism = name_match.group('name')
        if collect and name_match and not name_match.group('name') == organism:
            idx += 1
            ts = TabooSet(data, idx)
            connected = ts.connected
            logger.info("idx: %s | name: %s | taboo count: %s | %s", idx, organism, len(ts.taboo_strings),
                'not connected' if not connected else 'connected')
            organism = name_match.group('name')
            data = []
            if not connected:
                not_connected_count += 1
        elif rsite_match and collect:
            restriction_site = rsite_match.group('nucleotides')
            #XXX but lot of N's cause memory error, we could of course enumerate them on the fly
            if sum(1 for i in restriction_site if i == SPECIAL_CHARACTER_N) <= max_wildcard:
                data.append(restriction_site)
            collect = False

    if data:
        logger.debug("last one")
        idx += 1
        ts = TabooSet(data, idx)
        connected = ts.connected
        logger.info("idx: %s | name: %s | taboo count: %s | %s", idx, organism, len(ts.taboo_strings),
            'not' if not connected else '' + 'connected')

    logger.info("progress: %d", idx)
    logger.info("finished in %ss", time.time() - start_time)
    logger.info("not connected count: %s", not_connected_count)



if __name__ == "__main__":
    description = """
Determine connectivity for Hamming graphs based on taboo sets.

recommended usage:
cat input.txt | python connectivity.py
"""

    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter
        )
    parser.add_argument(
        '--organism-regex', default='<3>(?P<name>.*)',
        )
    parser.add_argument(
        '--enzyme-regex', default='^<1>.*[^I]II[AB]?$',
        )
    parser.add_argument(
        '--rsite-regex', default='<5>(?P<nucleotides>.*)',
        )
    parser.add_argument(
        '--max-wildcard', default=4, type=int,
        )
    parser.add_argument(
        '--log-level', default='DEBUG',
        )
    args = parser.parse_args()

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    ch.setLevel(args.log_level)
    logger.addHandler(ch)

    run(args.organism_regex, args.enzyme_regex, args.rsite_regex, args.max_wildcard)
