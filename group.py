from collections import defaultdict, Counter

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
        components = sorted(components, key=lambda x: len(x))
        component_sizes = [len(c) for c in components]
        if not cc == 2:
            print(f"Taboos: {taboos}  |  Components: {components}  |  Sizes: {component_sizes}")
            assert False, \
            f"Number of connected components is {cc}"
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

    def _generate_taboo_count_sums(self):
        for taboos, idx in self.orbit_map.items():
            self.index_map[idx].append(taboos)
        for orbit_idx, (c1, c2, t_count) in self.size_orbit_map.items():
            sum_taboo_count = t_count + c1
            self.taboo_count_sums[sum_taboo_count].append((1, orbit_idx, 0))
            sum_taboo_count = t_count + c2
            self.taboo_count_sums[sum_taboo_count].append((1, orbit_idx, 1))
        for (orbit_idx_1, (c1_1, c2_1, t_count_1)), (orbit_idx_2, (c1_2, c2_2, t_count_2)) \
            in combinations_with_replacement(self.size_orbit_map.items(), 2):
            sum_taboo_count = t_count_1 + t_count_2
            self.taboo_count_sums[sum_taboo_count].append((0, orbit_idx_1, orbit_idx_2))

    def extend(self):
        hc = self.__class__(self.n+1)
        self._generate_taboo_count_sums()
        for taboo_count in range(self.n+1, 2**(self.n) + 1):
            patterns = self.taboo_count_sums.get(taboo_count)
            if not patterns:
                continue
            for pattern in patterns: 
                if pattern[0] == 1:
                    orbit_idx = pattern[1]
                    comp_idx = pattern[2]
                    taboos_list = self.index_map[orbit_idx]
                    # cannot be connected
                    #connected_count = 0
                    extension_count = 0
                    for t_idx, taboos in enumerate(taboos_list, 1):
                        component = self.components[taboos][comp_idx]
                        next_taboos = [tuple(list(t) + [0]) for t in taboos]
                        next_taboos += [tuple(list(t) + [1]) for t in component]
                        next_taboos = frozenset(next_taboos)
                        if any(t.issubset(next_taboos) for t in hc.orbit_map):
                            extension_count += 1
                            continue
                        dc = hc.explore(next_taboos)
                        next_orbit_idx = hc.orbit_map[next_taboos]
                        hc.parents[next_orbit_idx] = (orbit_idx, comp_idx)
                    print(f"Pattern: {pattern} | Total: {t_idx}  |  Extension count: {extension_count}")
                if pattern[0] == 0:
                    orbit_idx_1 = pattern[1]
                    orbit_idx_2 = pattern[2]
                    taboos_list_1 = self.index_map[orbit_idx_1]
                    if orbit_idx_1 == orbit_idx_2:
                        products = combinations_with_replacement(taboos_list_1, 2)
                    else:
                        taboos_list_2 = self.index_map[orbit_idx_2]
                        products = direct_product([taboos_list_1, taboos_list_2])
                    connected_count = 0
                    extension_count = 0
                    for t_idx, (taboos_1, taboos_2) in enumerate(products, 1):
                        next_taboos = [tuple(list(t) + [0]) for t in taboos_1]
                        next_taboos += [tuple(list(t) + [1]) for t in taboos_2]
                        next_taboos = frozenset(next_taboos)
                        # prefilter 1: extension
                        if any(t.issubset(next_taboos) for t in hc.orbit_map):
                            extension_count += 1
                            continue
                        # prefilter 2: connected
                        if any(all(c1.intersection(c2) 
                            for c2 in self.components[taboos_2])
                            for c1 in self.components[taboos_1]):
                            connected_count += 1
                            continue
                        dc = hc.explore(next_taboos)
                        assert dc
                        next_orbit_idx = hc.orbit_map[next_taboos]
                        hc.parents[next_orbit_idx] = \
                            (orbit_idx_1, orbit_idx_2, int(taboos_1 == taboos_2))
                    print(f"Pattern: {pattern} | Total: {t_idx}  |  Connected count: {connected_count}  |  Extension count: {extension_count}")
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
