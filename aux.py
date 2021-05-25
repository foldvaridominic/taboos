import logging
import random
from collections import Counter, defaultdict
from functools import reduce
from itertools import combinations
import networkx as nx


logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
ch.setLevel(logging.INFO)
logger.addHandler(ch)


from constants import NUCLEOTIDE_CHARACTERS, ALPHABET_LENGTH
from tabooset import TabooSet as TS
from utils import *
from plot import DrawCube


nucleotides = set(NUCLEOTIDE_CHARACTERS)
allowed = []


def rec_func(length, letters=nucleotides):
    global allowed
    for r in range(1,len(letters) + 1):
        allowed.append(set(random.sample(letters, r)))
        if length > 1:
            yield from rec_func(length-1, letters)
        else:
            yield allowed
        allowed = allowed[:-1]


def gen_hamming_graph(taboos, length, letters=nucleotides):
        graph = nx.Graph()
        nodes = list(get_self_product(letters, length))
        edges = combinations(nodes, 2)
        edges = [e for e in edges if hamming_distance_1_for_strings(e)]
        graph.add_nodes_from(nodes)
        graph.add_edges_from(edges)
        graph.remove_nodes_from(taboos)
        cc = list(nx.algorithms.components.connected_components(graph))
        blocks = [len(c) for c in cc]
        #print(f"Connected components: {cc}")
        print(f"Size of connected components: {blocks}")
        print(f"Number of connected components: {len(blocks)}")
        return blocks


def create_projections(taboos, length, letters=nucleotides):
    cst = []
    alphabet = len(letters)
    for fixed in range(length):
        projection = [i for i in range(length) if i != fixed]
        for a in letters:
            projected_taboos = [tuple(t[i] for i in projection) 
                for t in taboos if t[fixed] == a]
            dummy_graph = nx.Graph()
            dummy_nodes = list(get_self_product(letters, length-1))
            dummy_graph.add_nodes_from(dummy_nodes)
            dummy_edges = combinations(dummy_nodes, 2)
            dummy_edges = [e for e in dummy_edges if hamming_distance_1_for_strings(e)]
            dummy_graph.add_edges_from(dummy_edges)
            dummy_graph.remove_nodes_from(projected_taboos)
            taboo_edges = [[e for e in direct_product([list(dummy_graph.nodes), [n]])
                      if hamming_distance_1_for_strings(e)]
                      for n in projected_taboos]
            cc = nx.algorithms.components.number_connected_components(dummy_graph)
            if cc > 1:
                minimal = True
                for te in taboo_edges:
                    dg = dummy_graph.copy()
                    dg.add_edges_from(te)
                    if nx.algorithms.components.number_connected_components(dg) > 1:
                            minimal = False
                            break
                if minimal:
                    cst.append("MC")
                else:
                    cst.append("DC")
            else:
                cst.append("C")
    cst = frozenset(Counter([frozenset(Counter([cst[k+i] for i in range(alphabet)]).items()) for k in range(0, len(cst), alphabet)]).items())
    return cst



def enum_graphs(length, alphabet=None):
    letters = set(range(alphabet)) if alphabet else nucleotides
    all_structures = defaultdict(list)
    all_projections = defaultdict(list)
    for allowed in rec_func(length, letters):
        not_allowed = [letters-pos for pos in allowed]
        taboos = []
        for idx in range(length):
            taboos += [s for s in direct_product(
                allowed[:idx] + [not_allowed[idx]] + allowed[idx+1:])]
        structure = tuple([len(a) for a in allowed])
        print(f"Taboo blocks: {structure}")
        blocks = gen_hamming_graph(taboos, length, letters)
        if len(blocks) > 1:
            block = frozenset(Counter(blocks).items())
            cst = create_projections(taboos, length, letters)
            all_structures[block].append(structure)
            all_projections[cst].append(structure)
    print(f"Number of different block structures: {len(all_structures)}")
    print(f"All different block structures")
    for block, count in all_structures.items():
        print(f"{block}: {count}")
    print(f"Number of different proj structures: {len(all_projections)}")
    print(f"All different proj structures")
    for cst, count in all_projections.items():
        print(f"{cst}: {count}")


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
        blocks = gen_hamming_graph(taboos, idx+1+jdx)
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


def inspect_dimension_increment_in_quotient_graph(characters, increase=1):
    node_groups, index_groups = enum_quotient_graph(characters)
    for n in node_groups:
        logger.info(n)
    for i in range(increase):
        node_groups, index_groups = overlaps(node_groups, index_groups)
        for n in node_groups:
            logger.info(n)


class TabooTree:

    def __init__(self, length, alphabet):
        nodes = list(get_self_product(range(alphabet), length))
        graph = nx.Graph()  
        graph.add_nodes_from(nodes)
        edges = combinations(nodes, 2)
        edges = [e for e in edges if hamming_distance_1_for_strings(e)]
        graph.add_edges_from(edges)
        #logger.info("Start: %s", list(graph.edges))
        self.graph = graph
        current_branch = [{n,} for n in nodes]
        self.current_branch = current_branch
        self.graph_nodes = nodes
        self.num_states = alphabet**length
        self.draw_cube = DrawCube(alphabet, length)
        self.n = length
        self.alphabet = range(alphabet)
        self.check_connected()

    def create_new_branch(self, length, skip):
        for c in combinations(self.graph_nodes, length):
            if not any(s <= set(c) for s in skip):
                yield set(c)

    def check_connected(self):
        new_branch = self.current_branch
        disconnected = []
        for i in range(self.num_states-1):
            # node in the TabooTree is a collection of nodes to be removed from the Hamming-graph
            component_sizes = defaultdict(int)
            cross_section_types = defaultdict(int)
            count = 0
            for idx, remove_nodes in enumerate(new_branch, 1):
                remove_edges = list(flatten([list(self.graph.edges(n)) for n in remove_nodes]))
                #logger.info("%s Remove nodes: %s", idx, remove_nodes)
                #logger.info("%s Remove edges: %s", idx, remove_edges)
                self.graph.remove_edges_from(remove_edges)
                self.graph.remove_nodes_from(remove_nodes)
                #logger.info("%s Remaining edges: %s", idx, list(self.graph.edges))
                cc = nx.algorithms.components.number_connected_components(self.graph)
                components = list(nx.algorithms.components.connected_components(self.graph))
                self.graph.add_edges_from(remove_edges)
                if cc > 1:
                    disconnected.append(remove_nodes)
                    component_size = frozenset([len(c) for c in components])
                    #logger.info("Taboo count: %s | idx: %s | components: %s", i+1, idx, component_size)
                    component_sizes[component_size] += 1
                    count += 1
                    #if self.num_states == 16:
                        #cst = self.draw_cube.create_fig_with_projections(remove_nodes, f'{i+1}_{count}')
                        # resize and regroup
                        #cst = frozenset(Counter([frozenset((cst[k],cst[k+1])) for k in range(0, len(cst), 2)]).items())
                        #cross_section_types[cst] += 1
                    cst = create_projections(remove_nodes, self.n, self.alphabet)
                    cross_section_types[cst] += 1
            logger.info("Taboo count %s finished: %s", i+1, count)
            logger.info(component_sizes.items())
            for cst, count in cross_section_types.items():
                logger.info("%s: %s", cst, count)
            new_branch = self.create_new_branch(i+2, disconnected)
        #logger.info("End: %s", list(self.graph.edges))
