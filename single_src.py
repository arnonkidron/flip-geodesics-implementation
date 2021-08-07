from path_shortener import PathShortener
from numeric_error_thresholds import REFLEX_ANGLE_THRESHOLD_FOR_SINGLE_SRC_GEODESICS
from queue import PriorityQueue
import numpy as np


class SingleSrcShortener:
    def __init__(self, triangulation):
        self.tri = triangulation
        self.src = None
        self.parent = None
        self.visited = None
        self.d = None
        self.q = None

        self.reflex_angle_threshold_for_edge_flip = REFLEX_ANGLE_THRESHOLD_FOR_SINGLE_SRC_GEODESICS
        self.reflex_angle_threshold_for_flip_out = REFLEX_ANGLE_THRESHOLD_FOR_SINGLE_SRC_GEODESICS

    def get_path(self):
        if self.parent is None:
            return [self.src]
        return [[v, self.parent[v]]
                for v in range(len(self.tri.V))
                if not np.isinf(self.d[v]) and v != self.src and self.visited[v] and self.parent[v] is not None]

    def get_frontier(self):
        return [[v, self.parent[v]]
                for v in range(len(self.tri.V))
                if not np.isinf(self.d[v]) and v != self.src and not self.visited[v]]

    def set_path(self, source_vertex_index):
        if isinstance(source_vertex_index, list):
            source_vertex_index = source_vertex_index[0]
        self.src = source_vertex_index

    def shorten(self, path):
        """
        apply MakeGeodesic on the path
        """
        shortener = PathShortener(self.tri)
        shortener.reflex_angle_threshold_for_edge_flip = self.reflex_angle_threshold_for_edge_flip
        shortener.reflex_angle_threshold_for_flip_out = self.reflex_angle_threshold_for_flip_out
        shortener.keep_list_of_flipped_edges = True
        shortener.set_path(path)
        shortener.make_geodesic()
        return shortener.length, shortener.get_path(), shortener.flipped_edges

    def make_geodesic(self, limit_iterations=None, length_threshold=None):
        for _ in self.make_geodesic_one_at_a_time(limit_iterations, length_threshold):
            pass

    def make_geodesic_one_at_a_time(self, limit_iterations=None, length_threshold=None):
        """
        Create the triangulation that includes all geodesic paths from the single source
        """
        num_v = len(self.tri.V)
        src = self.src
        self.d = np.ones(num_v) * np.inf
        self.visited = np.zeros(num_v, dtype=bool)
        self.parent = [None for _ in range(num_v)]

        self.d[src] = 0

        self.q = PriorityQueue(maxsize=num_v)
        self.q.put((self.d[src], src))

        while not self.q.empty():
            _, v = self.q.get()
            if self.visited[v]:
                continue
            self.visited[v] = True

            path_to_v = self.tri.find_shortest_path(self.src, v, self.parent)

            for e in self.tri.in_edges[v]:
                u = e.origin
                if self.visited[u]:
                    continue
                possible_dist, short_path, flipped_edges = self.shorten(path_to_v + [u])
                self.reconnect_vertices(flipped_edges)

                if self.d[u] > possible_dist:
                    self.d[u] = possible_dist
                    self.q.put((self.d[u], u))

                    for i in range(len(short_path) - 1):
                        self.parent[short_path[i+1]] = short_path[i]

            yield

    def reconnect_vertices(self, flipped_edges):
        for e in flipped_edges:
            u, v = e.origin, e.dst
            if self.parent[u] is not None and self.parent[u] == v:
                self.reconnect_vertex(u)
                continue
            v, u = u, v
            if self.parent[u] is not None and self.parent[u] == v:
                self.reconnect_vertex(u)

    def reconnect_vertex(self, u):
        self.d[u] = np.inf
        self.parent[u] = None
        for e in self.tri.in_edges[u]:
            v = e.origin
            possible_dist = self.d[v] + e.length
            if self.d[u] > possible_dist:
                self.d[u] = possible_dist
                self.q.put((self.d[u], u))
                self.parent[u] = v
