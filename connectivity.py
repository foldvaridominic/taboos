import sys
import argparse
import logging
import time
import random
import itertools
import re

import networkx as nx

from constants import *


def flatten(iterable):
    return list(itertools.chain.from_iterable(iterable))


class TabooSet:
    def __init__(self, taboos):
        logger.info("taboo init data: %s", taboos)
        first_last_letters = set()
        for t in taboos:
            first_last_letters.add(ALL_CHARACTERS.get(t[0]))
            first_last_letters.add(NUCLEOTIDE_COMPLEMENTS.get(t[-1]))
        first_last_letters = set(flatten([l for l in first_last_letters if l]))
        #logger.debug("first and last: %s", first_last_letters)
        self.taboo_init_counts = len(first_last_letters)
        self.taboo_inits = first_last_letters
        self.taboo_strings = taboos
        if self.taboo_init_counts >= ALPHABET_LENGTH:
            self.taboos = set(self.transform(taboos))
            #logger.debug("transformed taboos: %s", self.taboos)
            self.complement()
            self.minimize()
            self.taboo_strings = {''.join(t) for t in self.taboos}
            self.taboo_inits = {t[0] for t in self.taboos if t}
            self.taboo_init_counts = len(self.taboo_inits)
            self.max_taboo_length = max(map(len, self.taboo_strings))
            self.min_taboo_length = min(map(len, self.taboo_strings))
            self.nodes = dict()
            self.graph = nx.Graph()

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
        sorted_taboos = set(sorted(
            [t[1:] for t in self.taboo_strings
            if all(tt not in t[1:] for tt in self.taboo_strings)],
            key=len, reverse=True))
        pairs = list(itertools.combinations_with_replacement(sorted_taboos, 2))
        #logger.debug("Taboo pairs: %s", pairs)
        pairs_filtered = [p for p in pairs if sum(1 for t1,t2 in zip(*p) if t1 != t2) <= 1]
        #logger.debug("Taboo pairs filtered: %s", pairs_filtered)
        logger.debug("Taboo pairs filtered count: %s", len(pairs_filtered))
        for pair in pairs_filtered:
            for idx, letter in enumerate(NUCLEOTIDE_CHARACTERS, 1):
                first = letter + pair[0]
                second = letter + pair[1]
                if all(t not in first and t not in second for t in self.taboo_strings):
                    break
                else:
                    #logger.debug("%s or %s in taboo set", first, second)
                    pass
                if idx == ALPHABET_LENGTH:
                    logger.debug("pair %s cannot be left-1 synchronized", pair)
                    return False
        return True

    def connected_3(self):
        self.gen_nodes_with_length(self.max_taboo_length)
        self.gen_suffix_classifiers()
        #logger.debug(self.nodes)
        logger.debug("Taboo strings: %s", self.taboo_strings)
        logger.debug("LSC: %s", self.long_suffix_classifiers)
        logger.debug("SSC: %s", self.short_suffix_classifiers)
        M_length_nodes = self.nodes[self.max_taboo_length]
        for r in self.long_suffix_classifiers:
            logger.info("Checking long suffix classifier: %s", r)
            if not self.gen_suffix_graph(r, len(r)+self.max_taboo_length, M_length_nodes,
                self.max_taboo_length, quotient=True, partition_length=1):
                return False
        for p in self.short_suffix_classifiers:
            p_length = len(p)
            partition_length = self.max_taboo_length - p_length
            logger.info("Checking short suffix classifier: %s", p)
            if not self.gen_suffix_graph(p, 2*self.max_taboo_length-1, M_length_nodes,
                    self.max_taboo_length, quotient=True, partition_length=partition_length):
                return False
            for length in range(p_length+2, self.max_taboo_length):
                initial_nodes = self.nodes.get(length)
                if not initial_nodes:
                    continue
                if not self.gen_suffix_graph(p, length, initial_nodes, self.max_taboo_length):
                    return False
        return True

    def get_initial_nodes(self, starting_point):
        initial_nodes = [''.join(t) for t in itertools.product(NUCLEOTIDE_CHARACTERS, repeat=starting_point)]
        return initial_nodes

    def gen_nodes_with_length(self, length, initial_nodes=None, save=True):
        starting_point = len(initial_nodes[0]) if initial_nodes else self.min_taboo_length - 1
        current_nodes = initial_nodes or self.get_initial_nodes(starting_point)
        while starting_point < length:
            current_nodes = self.extend_nodes(current_nodes)
            starting_point += 1
            if save:
                self.nodes[starting_point] = current_nodes
        return current_nodes

    def extend_nodes(self, current_nodes):
        ret = []
        M_length = len(current_nodes[0])
        for node in current_nodes:
            ext = [NUCLEOTIDE_CHARACTERS[idx] + node
            for idx in range(ALPHABET_LENGTH)
            if all(t not in NUCLEOTIDE_CHARACTERS[idx] + node
            for t in self.taboo_strings)]
            if M_length and not ext: 
                logger.info("Not a left proper taboo set")
            ret += ext
        logger.debug("extending nodes to length %s", len(current_nodes[0])+1)
        return ret

    def gen_suffix_classifiers(self):
        M_length_nodes = self.nodes[self.max_taboo_length]
        taboo_suffixes = set([t[idx::] for t in self.taboo_strings for idx in range(len(t),0,-1)])
        M_suffixes = set([tf[idx::] for tf in M_length_nodes for idx in range(len(tf),0,-1) ])
        scs = taboo_suffixes.intersection(M_suffixes)
        lscs = set()
        for node in M_length_nodes:
            re_node = re.compile(node)
            lsc = list(zip(*sorted(set([(sc, re_node.match(sc).end())
            if re_node.match(sc) else ('', 0) for sc in scs]),
            key=lambda x: x[1], reverse=True)))[0][0]
            lscs.add(lsc)
        sscs = scs - lscs
        self.long_suffix_classifiers = lscs
        self.short_suffix_classifiers = sscs
        return lscs, sscs

    def gen_suffix_graph(self, suffix, node_length, initial_nodes,
            initial_length, quotient=False, partition_length=None):
        self.graph.clear()
        if initial_length >= len(suffix):
            initial_nodes = [n for n in initial_nodes if n.endswith(suffix)]
        if initial_length < node_length:
            current_nodes = self.gen_nodes_with_length(node_length, 
                initial_nodes=initial_nodes, save=False)
        else:
            current_nodes = initial_nodes
        if quotient:
            edges = []
            nodes = []
            partitions = self.nodes.get(partition_length) or \
                [''.join(p) for p in itertools.product(NUCLEOTIDE_CHARACTERS, repeat=partition_length)]
            partitions = [p for p in partitions if any(n.endswith(p+suffix) for n in current_nodes)]
            remaining_length = node_length - len(suffix) - partition_length
            logger.debug("Partition size: %s", len(partitions))
            for i1, p1 in enumerate(partitions, 1):
                nodes.append(p1)
                for p2 in partitions[i1:]:
                    if sum(1 for t1,t2 in zip(p1,p2) if t1 != t2) > 1:
                        continue
                    prefixes = {n[:remaining_length] for n in current_nodes}
                    if any(prefix+p1+suffix in current_nodes and
                        prefix+p2+suffix in current_nodes for prefix in prefixes):
                        edges.append((p1, p2))
        else:
            edges = itertools.combinations(current_nodes, 2)
            edges = [p for p in edges if sum(1 for t1,t2 in zip(*p) if t1 != t2) <= 1]
            nodes = current_nodes
        self.graph.add_nodes_from(nodes)
        self.graph.add_edges_from(edges)
        cc = nx.algorithms.components.number_connected_components(self.graph)
        logger.info(
            "Connected components of %s suffix graph %s with length %s: %s",
            'quotient ' + str(partition_length) if quotient else '', suffix, node_length, cc
            )
        return cc == 1

    @property
    def connected(self):
        if self.connected_1:
            logger.debug("Connected 1: %s", self.taboo_inits)
            return True
        if self.connected_2:
            logger.debug("Connected 2")
            return True
        return self.connected_3()


def run(organism_regex, enzyme_regex, rsite_regex, max_wildcard):

    idx=0
    start_time = time.time()
    re_organism = re.compile(organism_regex)
    re_enzyme = re.compile(enzyme_regex)
    re_rsite = re.compile(rsite_regex)
    organism = None
    prev_org = None
    data = []
    collect = False
    not_connected_count = 0


    for line in sys.stdin:
        name_match = re_organism.search(line)
        enzyme_match = re_enzyme.search(line)
        rsite_match = re_rsite.search(line)

        if enzyme_match:
            collect = True
        if collect and name_match:
            prev_org = organism if organism is not None else name_match.group('name')
            organism = name_match.group('name')
        if collect and data and prev_org != organism:
            idx += 1
            logger.info("%s name: %s", idx, prev_org)
            ts = TabooSet(data)
            connected = ts.connected
            logger.info("taboo count: %s | %s", len(ts.taboo_strings),
                'not connected' if not connected else 'connected')
            data = []
            prev_org = organism
            if not connected:
                not_connected_count += 1
        if collect and rsite_match:
            restriction_site = rsite_match.group('nucleotides')
            if sum(1 for i in restriction_site if i == SPECIAL_CHARACTER_N) <= max_wildcard:
                current_data = restriction_site
            else:
                current_data = 'too many wildcard characters'
            data.append(current_data)
            collect = False
            prev_org = organism

    if data:
        idx += 1
        logger.info("%s name: %s", idx, prev_org)
        logger.info("last one")
        ts = TabooSet(data)
        connected = ts.connected
        logger.info("taboo count: %s | %s", len(ts.taboo_strings),
            'not connected' if not connected else 'connected')

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
        '--organism-regex', default='^OS\s+(?P<name>[^\r]*)',
        )
    parser.add_argument(
        '--enzyme-regex', default='^ET\s+RM?2',
        )
    parser.add_argument(
        '--rsite-regex', default='^RS\s+(?P<nucleotides>[^,]*),',
        )
    parser.add_argument(
        '--max-wildcard', default=8, type=int,
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
