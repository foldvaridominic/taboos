import time
from collections import defaultdict, Counter

import networkx as nx
from sympy.combinatorics.named_groups import SymmetricGroup as SG

from utils import *


class HGSymmetryGroup:

    def __init__(self, n, q):
        '''
        Represent the symmetry group of n-dimensional Hamming-graph
        as the direct product of n S_q and S_n. 
        '''
        self.q_permutations = SG(q)
        self.n_permutations = SG(n)
        self.q = q
        self.n = n
        self.mirrors = Mirrors(list(self.q_permutations.generate()), n)
        self.alphabet = range(q)
        self.dimension = range(n)

    def orbit(self, points):
        '''
        :points: frozenset of tuples of 0s and 1s of length dimension
        return: the orbit of the points
        '''
        counts = set()
        for mirrors in self.mirrors.generate():
            q_action = [[mirrors[i](self.alphabet)[point[i]] for i in self.dimension]
                for point in points]
            for permutation in self.n_permutations.generate():
                ret = frozenset(tuple(permutation(qa)) for qa in q_action)
                if ret in counts:
                    continue
                counts.add(ret)
                yield ret


class HammingGraph:

    def __init__(self, n, q):
        self.n = n
        self.q = q
        self.num_states = q**n
        graph = nx.Graph()
        nodes = list(get_self_product(range(q), n))
        graph.add_nodes_from(nodes)
        edges = combinations(nodes, 2)
        edges = [e for e in edges if hamming_distance_1_for_strings(e)]
        graph.add_edges_from(edges)
        self.graph = graph
        self.group = HGSymmetryGroup(n,q)
        self.orbit_map = defaultdict(int)
        self.size_orbit_map = {}
        self.components = {}
        self.parents = {}
        self.taboo_count_sums = defaultdict(list)
        self.index_map = defaultdict(list)


    def disconnect(self, taboos):
        remove_nodes = list(taboos)
        remove_edges = list(flatten([list(self.graph.edges(n)) for n in taboos]))
        #print(f"Remove nodes: {remove_nodes}")
        #print(f"Remove edges: {remove_edges}")
        self.graph.remove_edges_from(remove_edges)
        self.graph.remove_nodes_from(remove_nodes)
        cc = nx.algorithms.components.number_connected_components(self.graph)
        components = nx.algorithms.components.connected_components(self.graph)
        if cc == 1:
            self.graph.add_edges_from(remove_edges)
            return False
        components = list(map(lambda x: frozenset(x), sorted(components, key=lambda x: len(x))))
        component_sizes = [len(c) for c in components]
        #if not cc == 2:
        #    print(f"Taboos: {taboos}  |  Components: {components}  |  Sizes: {component_sizes}")
        #    assert False, \
        #    f"Number of connected components is {cc}"
        self.graph.add_edges_from(remove_edges)
        return component_sizes, components


    def explore(self, taboos=[(0,0),(1,1)]):
        orbit_idx = max(self.orbit_map.values()) + 1 if self.orbit_map else 1
        orbit = self.group.orbit(taboos)
        for idx, taboos in enumerate(orbit, 1):
            if taboos in self.orbit_map:
                print(f"Orbit {self.orbit_map[taboos]} already generated")
                return False
            disconnect = self.disconnect(taboos)
            taboo_count = len(taboos)
            if idx == 1 and disconnect:
                component_sizes, components = disconnect 
                self.size_orbit_map[orbit_idx] = tuple(component_sizes + [taboo_count])
            elif disconnect:
                sizes, components = disconnect
                assert component_sizes == sizes, \
                    f"Component sizes changed from {component_sizes} to {sizes}"
            else:
                return False
            self.orbit_map[taboos] = orbit_idx
            self.components[taboos] = components
        print(f"{orbit_idx}. Orbit size: {idx} | Components: {component_sizes} | Taboo count: {taboo_count}")
        return True

    def create_new_branch(self, length):
        for c in combinations(self.graph.nodes, length):
            yield frozenset(c)

    def extend(self):
        start_time = time.time()
        hc = self.__class__(self.n+1, self.q)
        for taboos_1 in self.orbit_map:
            orbit_idx_1 = self.orbit_map[taboos_1]
            components = self.components[taboos_1]
            extension_count = 0
            connected_count = 0
            new_orbit_count = 0
            for c in combinations_with_replacement(range(1,self.num_states+1),self.n):
                branches = [list(self.create_new_branch(i)) for i in c]
                all_taboos = direct_product(branches)
                for taboos in all_taboos:
                    next_taboos = [tuple(list(t) + [0]) for t in taboos_1]
                    for i in range(self.n):
                        next_taboos += [tuple(list(t) + [i+1]) for t in taboos[i]]
                    next_taboos = frozenset(next_taboos)
                    # prefilter 1: extension
                    if any(t.issubset(next_taboos) for t in hc.orbit_map):
                        extension_count += 1
                        continue
                    dc = hc.explore(next_taboos)
                    # prefilter 2: connected
                    if not dc:
                        connected_count +=1 
                        continue
                    next_orbit_idx = hc.orbit_map[next_taboos]
                    new_orbit_count += 1
                    orbit_idx_i = [(self.orbit_map.get(taboos[i], "U"), \
                        f"C{components.index(taboos[i])+1}" if taboos[i] in components \
                        else "S" if (taboos[i] == taboos_1) else "")
                        for i in range(self.n)]
                    hc.parents[next_orbit_idx] = (orbit_idx_1, orbit_idx_i)
                    print(f"Orbit: {next_orbit_idx} | {[orbit_idx_1] + orbit_idx_i}")
            print(f"Pattern: {orbit_idx_1}: {taboos_1} | New orbit: {new_orbit_count} | Connected count: {connected_count} | Extension count: {extension_count}")
        finished = time.time() - start_time
        print(f"Finished in: {finished} s")
        return hc

    def create_projections(self): 
        orbit_projection_map = {}
        coord_indices = list(range(self.n))
        cross_section_map = defaultdict(set)
        for idx, (taboos, orbit_idx) in enumerate(self.orbit_map.items(), 1):
            if not idx % 1000: 
                print(f"Progress: {idx}")
            cst = []
            for fixed in coord_indices:
                projection = [i for i in range(self.n) if i != fixed]
                for a in range(2):
                    projected_taboos = [tuple(t[i] for i in projection) for t in taboos if t[fixed] == a]
                    dummy_graph = nx.Graph()
                    dummy_nodes = list(get_self_product(range(2), self.n-1))
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
            cst = frozenset(Counter([frozenset((cst[k],cst[k+1])) for k in range(0, len(cst), 2)]).items())
            cross_section_map[orbit_idx].add(cst)
        for orbit_idx, cross_section_types in cross_section_map.items():
            print(f"Orbit: {orbit_idx} | {cross_section_types}")
