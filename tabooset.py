import logging

import networkx as nx

from constants import *
from utils import *


logger = logging.getLogger()


class TabooSet:
    '''
    Check if Hamming-graph (n,4) is connected for all n with a certain taboo set.
    '''

    def __init__(self, taboos, complement=True, ignore_first_check=False):
        '''
        :param taboos: list of nucleotide seqs, wild card characters allowed
            Example: taboos = [GCTNS, AKYAGTB, TCCC]
        :param complement: boolean to consider complementary strands as well
        :param ignore_first_check: boolean to ignore first connectivity check
        '''

        logger.info("taboo init data: %s", taboos)

        # check if the whole alphabet is used as the first character
        first_last_letters = set()
        for t in taboos:
            first_last_letters.add(ALL_CHARACTERS.get(t[0]))
            first_last_letters.add(NUCLEOTIDE_COMPLEMENTS.get(t[-1]))
        first_last_letters = set(flatten([l for l in first_last_letters if l]))
        #logger.debug("first and last: %s", first_last_letters)

        self.taboo_init_counts = len(first_last_letters)
        self.taboo_strings = taboos
        if not self.connected_1 or ignore_first_check:
            taboo_tuples = self.transform(taboos)
            taboo_regexes = self.transform(taboos, regex=True)
            #logger.debug("transformed taboos: %s", taboo_regexes)
            complement_tuples = self.complement(taboo_tuples) if complement else set()
            complement_regexes = self.complement(taboo_regexes) if complement else set()
            #logger.debug("complement taboos: %s", complement_regexes)
            self.taboos = taboo_tuples | complement_tuples
            regexes = taboo_regexes | complement_regexes
            self.taboo_strings = {to_string(r) for r in regexes}
            logger.debug("all taboos: %s", self.taboo_strings)
            self.minimize()
            self.max_taboo_length = max(map(len, self.taboos)) if self.taboos else 1
            self.min_taboo_length = min(map(len, self.taboos)) if self.taboos else 1
            self._nodes = dict()
            self.graph = nx.Graph()

    def transform(self, taboos, regex=False):
        '''
        Parse the raw input nucleotide sequences.
        '''
        ret = set()
        for t in taboos:
            new = self.generate_chars(t, regex=regex)
            if new:
                ret.add(new)
        return ret

    def complement(self, taboos):
        ret = set()
        for t in taboos:
            ret.add(self.generate_complement(t))
        return ret

    def minimize(self):
        #TODO
        pass

    def generate_chars(self, taboo_string, regex=False):
        '''
        Parse a taboo sequence. Return a tuple where each element is
            regex=False: a tuple made from either a single valid nucleotide or 
                         a set nucleotides in case it is a wildcard character
            regex=True: a single valid nucleotide or a regular expression.
        When a character cannot be parsed ignore it.
        '''
        if regex:
            chars =  [l if l in NUCLEOTIDE_CHARACTERS 
                else SPECIAL_CHARACTERS_REGEX.get(l)
                for l in taboo_string]
        else:
            chars = [(l,) if l in NUCLEOTIDE_CHARACTERS 
                else SPECIAL_CHARACTERS.get(l)
                for l in taboo_string]

        # ignore all other characters for now
        chars_filtered = tuple(v for v in chars if v is not None)
        diff = len(chars) - len(chars_filtered)
        if diff:
            logger.debug("Unknown characters in %s", taboo_string)
        return chars_filtered

    def generate_complement(self, taboo_string):
        complement = tuple(NUCLEOTIDE_COMPLEMENTS[l] for l in taboo_string[::-1])
        return complement

    def taboo_free(self, instance):
        if isinstance(instance, tuple):
            func = taboo_free_for_tuples
        elif isinstance(instance, str):
            func = taboo_free_for_strings
        else:
            logger.error("Expected string or tuple, got %s", type(isinstance))
            raise TypeError
        return func(instance, self.taboo_strings)

    @property
    def connected(self):
        if self.connected_1:
            logger.debug("Connected 1: %s", self.taboo_init_counts)
            return True
        if self.connected_2:
            logger.debug("Connected 2")
            return True
        return self.connected_3()

    @property
    def connected_1(self):
        return self.taboo_init_counts < ALPHABET_LENGTH

    @property
    def connected_2(self):
        # Cut the first character from each taboo string
        sorted_taboos = set(sorted(
            [t[1:] for t in self.taboos
            if self.taboo_free(t[1:])],
            key=len, reverse=True))
        pairs = combinations_with_replacement(sorted_taboos, 2)
        #logger.debug("Taboo pairs: %s", pairs)
        pairs_filtered = [p for p in pairs if hamming_distance_1_for_tuples(p)]
        #logger.debug("Taboo pairs filtered: %s", pairs_filtered)
        logger.debug("Taboo pairs filtered count: %s", len(pairs_filtered))
        no_left_sync = []
        for p1, p2 in pairs_filtered:
            for product1 in direct_product(p1):
                product1_str = to_string(product1)
                if not self.taboo_free(product1_str):
                    continue
                for product2 in direct_product(p2):
                    product2_str = to_string(product2)
                    if not self.taboo_free(product2_str):
                        continue
                    if not hamming_distance_1_for_strings((product1_str, product2_str)):
                        continue
                    for idx, letter in enumerate(NUCLEOTIDE_CHARACTERS, 1):
                        first = letter + product1_str
                        second = letter + product2_str
                        if self.taboo_free(first) and self.taboo_free(second):
                            break
                        else:
                            #logger.debug("%s or %s in taboo set", first, second)
                            pass
                        if idx == ALPHABET_LENGTH:
                            logger.debug("pair %s and %s cannot be left-1 synchronized",
                                product1_str, product2_str)
                            no_left_sync.append((product1_str, product2_str))
        if no_left_sync:
            return False
        return True

    def connected_3(self):
        self.gen_nodes_with_length()
        self.gen_suffix_classifiers()
        #logger.debug(self._nodes)
        logger.debug("Taboo strings: %s", self.taboo_strings)
        logger.debug("LSC: %s", self.long_suffix_classifiers)
        logger.debug("SSC: %s", self.short_suffix_classifiers)
        for r in self.long_suffix_classifiers:
            logger.info("Checking long suffix classifier: %s", r)
            if not self.gen_suffix_graph(r, len(r)+self.max_taboo_length,
                quotient=True, partition_length=1):
                return False
        for p in self.short_suffix_classifiers:
            p_length = len(p)
            partition_length = max(self.min_lsc - p_length, 1)
            break_ = False
            while partition_length <= self.max_taboo_length - p_length:
                prefixes = self.get_nodes(partition_length)
                prefixes_filtered = [pf for pf in prefixes 
                    if self.taboo_free(pf+p)]
                # all possible prefix + ssc concat have to start with an lsc
                if all(any((pf+p).startswith(l)
                    for l in self.long_suffix_classifiers)
                    for pf in prefixes_filtered):
                    break
                partition_length += 1
            logger.info("Checking short suffix classifier %s with partition length %s",
                    p, partition_length)
            n = self.max_taboo_length - 1 + partition_length + p_length
            if not self.gen_suffix_graph(p, n, quotient=True,
                    partition_length=partition_length):
                return False
            for length in range(p_length+2, p_length + partition_length):
                initial_nodes = self._nodes.get(length)
                if not initial_nodes:
                    continue
                if not self.gen_suffix_graph(p, length, nodes=initial_nodes):
                    return False
        return True

    def get_initial_nodes(self, starting_point):
        return self.get_nodes(starting_point)

    def gen_nodes_with_length(self, length=None):
        '''
        :param length: integer
        :return: list of strings constructed as all taboo avoiding
                 character combinations up to length
        '''
        starting_point = self.min_taboo_length - 1
        length = length or self.max_taboo_length + 1
        current_nodes = self.get_initial_nodes(starting_point)
        while starting_point < length:
            current_nodes = self.extend_nodes(current_nodes)
            starting_point += 1
            logger.debug("extending nodes to length %s", starting_point)
            # cache nodes
            self._nodes[starting_point] = current_nodes
        return current_nodes

    def extend_nodes(self, current_nodes):
        '''
        :param current_nodes: list of strings
        :return: list of strings constructed as taboo avoiding left extension
                 of current nodes by all nucleotides
        '''
        ret = []
        try:
            M_length = len(current_nodes[0]) == self.max_taboo_length
        except (TypeError, IndexError):
            # initial nodes is a generator
            M_length = False
        for node in current_nodes:
            ext = [NUCLEOTIDE_CHARACTERS[idx] + node
            for idx in range(ALPHABET_LENGTH)
            if self.taboo_free(NUCLEOTIDE_CHARACTERS[idx] + node)]
            if M_length and not ext:
                logger.info("Not a left proper taboo set: %s", node)
            ret += ext
        return ret

    def gen_suffix_classifiers(self):
        M_length_nodes = self.get_nodes(self.max_taboo_length)
        # XXX memory hotspot
        taboo_suffixes = set([to_string(p[idx::]) for t in self.taboos 
            for p in direct_product(t) for idx in range(len(p),0,-1)])
        M_suffixes = set([tf[idx::] for tf in M_length_nodes for idx in range(len(tf),0,-1)])
        # suffix classifiers have to be suffixes of M-length taboo free strings
        scs = taboo_suffixes.intersection(M_suffixes)
        lscs = set()
        for node in M_length_nodes:
            lsc = list(zip(*sorted(set([(sc, match_from_start(sc, node).end())
            if match_from_start(sc, node) else ('', 0) for sc in scs]),
            key=lambda x: x[1], reverse=True)))[0][0]
            lscs.add(lsc)
        sscs = scs - lscs
        self.long_suffix_classifiers = lscs
        self.short_suffix_classifiers = sscs
        self.min_lsc = min(len(l) for l in lscs)
        return lscs, sscs

    def get_nodes(self, length=None): 
        '''
        :param length: integer
        :return: list or generator of taboo free strings up to length
        '''
        length = length or self.max_taboo_length - 1
        ret = self._nodes.get(length) or \
            (to_string(p) for p in get_self_product(NUCLEOTIDE_CHARACTERS, length))
        return ret

    def gen_suffix_graph(self, suffix, node_length, nodes=None,
            quotient=False, partition_length=None):
        '''
        Generate suffix graph.
        :param suffix: string
        :param node_length: integer
        :param nodes: list of strings
        :param quotient: boolean to generate quotient graph
        :param partition_length: integer, length of strings in quotient graph nodes
        :return: boolean whether the graph is connected or not
        '''
        self.graph.clear()
        if quotient:
            edges = []
            nodes = []
            # when node is represented as prefix + partition + suffix
            # both LSC and SSC have M-1 length prefixes
            partitions = list(self.get_nodes(partition_length))
            prefixes = self.get_nodes()
            for i1, p1 in enumerate(partitions, 1):
                p1_suffix = p1 + suffix
                if not any(self.taboo_free(prefix+p1_suffix) for prefix in prefixes):
                    # reset generator 
                    prefixes = self.get_nodes()
                    continue
                nodes.append(p1)
                for p2 in partitions[i1:]:
                    p2_suffix = p2 + suffix
                    if not hamming_distance_1_for_strings((p1,p2)):
                        continue
                    if any(self.taboo_free(prefix+p1_suffix) and
                        self.taboo_free(prefix+p2_suffix) for prefix in prefixes):
                        edges.append((p1, p2))
                    # reset generator 
                    prefixes = self.get_nodes()
        else:
            nodes = [n for n in nodes if n.endswith(suffix)]
            edges = combinations(nodes, 2)
            edges = [p for p in edges if hamming_distance_1_for_strings(p)]
        self.graph.add_nodes_from(nodes)
        self.graph.add_edges_from(edges)
        cc = nx.algorithms.components.number_connected_components(self.graph)
        logger.info(
            "Connected components of %s suffix graph %s with length %s: %s",
            'quotient ' + str(partition_length) if quotient else '', 
            suffix, node_length, cc
            )
        if cc > 1:
            logger.info(list(nx.algorithms.components.connected_components(self.graph)))
        return cc == 1
