import matplotlib.pyplot as plt
import numpy as np
import networkx as nx


from utils import *


class DrawCube:
    
    def __init__(self, alphabet, length):
        self.alphabet = range(alphabet)
        self.create_cube()
        self.length = length
        self.coord_indices = set(range(length))
        self.node_labels = {n: ''.join(map(str, n)) for n in list(self.graph)}

    def create_cube(self, x_offset=0.1, y_offset=3):
        g = nx.hypercube_graph(3)
        pos = nx.spring_layout(g)
        coordinates = [(0,1), (1,0), (-1,0), (0,-1)]
        coordinates = list(flatten([((c[0],c[1]),(c[0]+x_offset,c[1]+y_offset)) for c in coordinates]))
        for node, point in zip(list(g), coordinates):
            pos[node] = point
        self.graph = g
        self.pos = pos

    def create_node_colors(self, nodes):
        return ['black' if node in  nodes else 'white' for node in list(self.graph)]

    def create_projections(self, nodes, d=3):
        projections = combinations(self.coord_indices, d)
        for projection in projections:
            r_index = (self.coord_indices - set(projection)).pop()
            proj = ''.join(map(str, sorted(projection)))
            for a in self.alphabet:
                title = proj + '_' + str(a)
                yield [tuple(n[i] for i in projection) for n in nodes if n[r_index] == a], title

    def create_fig(self, nodes=[], filename='cube', path='./figs/'):
        fig, axes = plt.subplots(nrows=4, ncols=2)
        ax = axes.flatten()
        for idx, (group, title) in enumerate(nodes):
            node_colors = self.create_node_colors(group)
            nx.draw(self.graph, self.pos, ax=ax[idx], node_color=node_colors,
                labels=self.node_labels, font_color='yellow')
            ax[idx].set_title(title)
        path = path + filename
        fig.tight_layout()
        fig.savefig(path, dpi=fig.dpi)
        plt.close(fig)

    def create_fig_with_projections(self, nodes=[], filename='cube'):
        data =  [p for p in self.create_projections(nodes)]
        self.create_fig(data, filename)
