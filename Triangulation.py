from edge import Edge
import numpy as np
from copy import deepcopy

class Triangulation:
    def __init__(self, V, F):
        self.V = V
        self.edges = None
        self.init_edges(F)

    def init_edges(self, F):
        # initialize a list for every vertex
        self.edges = [[] for v in range(len(self.V))]

        # fill it
        for f in F:
            for i in range(3):
                u = f[i]
                v = f[(i+1) % 3]
                if u > v:
                    continue
                self.add_edge(u, v)
        return

    def add_edge(self, u, v):
        edge = Edge(u, v)
        self.edges[u].append(edge)
        self.edges[v].append(edge)
        return edge


    def remove_edge(self, u, v):
        for e in self.edges[u]:
            if e.endpoint1 == v or e.endpoint2 == v:
                self.edges[u].remove(e)

    def demo_flip(self):
        self.remove_edge(0,2)
        self.remove_edge(2,0)
        e = self.add_edge(1,3)

        e.midpoints.append([-7, 7, 7])


    def all_edges(self):
        for u in range(len(self.edges)):
            for e in self.edges[u]:
                v = e.endpoint1 if u == e.endpoint2 else e.endpoint2
                if u < v:
                    yield e

    def get_polyline(self):
        poly_vertices = deepcopy(self.V)
        poly_edges = []
        for e in self.all_edges():
            e_V, e_E = e.get_polyline(index_difference=len(poly_vertices))
            if e_V:
                poly_vertices = np.vstack([poly_vertices, e_V])
            poly_edges.extend(e_E)

        return poly_vertices, np.array(poly_edges)
