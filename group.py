from collections import defaultdict

import networkx as nx
from sympy.combinatorics.named_groups import SymmetricGroup as SG

from utils import *


class HyperOctahedralGroup:

    def __init__(self, n):
        '''
        Represent the n-dimensional hyperoctahedral group as a direct product of
        n perpendicular mirorrs and the symmetric group S_n. 
        '''
        self.mirrors = Mirrors(range(2), n)
        self.pos_permutations = SG(n)

    def orbit(self, points):
        '''
        :points: frozenset of tuples of 0s and 1s of length dimension
        return: the orbit of the points
        '''
        counts = set()
        for mirror in self.mirrors.generate():
            mirrored = [[i^m for i,m in zip(point, mirror)]
                for point in points]
            for permutation in self.pos_permutations.generate():
                ret = frozenset(tuple(permutation(m)) for m in mirrored)
                if ret in counts:
                    continue
                counts.add(ret)
                yield ret


class HyperCubeGraph:

    def __init__(self, n):
        self.n = n
        graph = nx.Graph()
        nodes = list(get_self_product(range(2), n))
        graph.add_nodes_from(nodes)
        edges = combinations(nodes, 2)
        edges = [e for e in edges if hamming_distance_1_for_strings(e)]
        graph.add_edges_from(edges)
        self.graph = graph
        self.group = HyperOctahedralGroup(n)
        self.orbit_map = defaultdict(int)
        self.size_orbit_map = {}
        self.components = {}
        self.parents = {}


    def disconnect(self, taboos):
        remove_nodes = list(taboos)
        remove_edges = list(flatten([list(self.graph.edges(n)) for n in taboos]))
        #print(f"Remove nodes: {remove_nodes}")
        #print(f"Remove edges: {remove_edges}")
        self.graph.remove_edges_from(remove_edges)
        self.graph.remove_nodes_from(remove_nodes)
        cc = nx.algorithms.components.number_connected_components(self.graph)
        if cc == 1:
            return False
        components = list(nx.algorithms.components.connected_components(self.graph))
        component_sizes = tuple(sorted(len(c) for c in components))
        self.graph.add_edges_from(remove_edges)
        assert cc == 2, f"Number of connected components is {cc}"
        #print(f"Taboos: {taboos}  |  Components: {cc}  |  Sizes: {component_sizes}")
        return component_sizes, components


    def explore(self, taboos=[(0,0),(1,1)]):
        orbit_idx = max(self.orbit_map.values()) + 1 if self.orbit_map else 1
        orbit = self.group.orbit(taboos)
        for idx, taboos in enumerate(orbit, 1):
            if taboos in self.orbit_map:
                print(f"Orbit {self.orbit_map[taboos]} already generated")
                return False
            disconnect = self.disconnect(taboos)
            if idx == 1 and disconnect:
                component_sizes, components = disconnect 
                self.size_orbit_map[orbit_idx] = component_sizes
            elif disconnect:
                sizes, components = disconnect
                assert component_sizes == sizes, \
                    f"Component sizes changed from {component_sizes} to {sizes}"
            else:
                return False
            self.orbit_map[taboos] = orbit_idx
            self.components[taboos] = components
        print(f"{orbit_idx}. Orbit size: {idx} | Components: {component_sizes}")
        return True

    def extend(self):
        hc = self.__class__(self.n+1)
        index_map = defaultdict(list)
        for taboos, idx in self.orbit_map.items():
            index_map[idx].append(taboos)
        for orbit_idx, taboos_list in index_map.items():
            for taboos in taboos_list:
                for c_idx, c in enumerate(self.components[taboos], 1):
                    next_taboos = [tuple(list(t) + [0]) for t in taboos]
                    next_taboos += [tuple(list(t) + [1]) for t in c]
                    next_taboos = frozenset(next_taboos)
                    dc = hc.explore(next_taboos)
                    next_orbit_idx = hc.orbit_map[next_taboos]
                    hc.parents[next_orbit_idx] = (orbit_idx, c_idx)
        for t1, t2 in combinations_with_replacement(self.orbit_map.keys(), 2):
            next_taboos = [tuple(list(t) + [0]) for t in t1]
            next_taboos += [tuple(list(t) + [1]) for t in t2]
            next_taboos = frozenset(next_taboos)
            dc = hc.explore(next_taboos)
            if not dc:
                continue
            next_orbit_idx = hc.orbit_map[next_taboos]
            t1_idx = self.orbit_map[t1]
            t2_idx = self.orbit_map[t2]
            hc.parents[next_orbit_idx] = (t1_idx, t2_idx, int(t1 == t2))
        return hc
