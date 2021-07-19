import numpy as np


class FlipEdgeNetwork:
    def __init__(self, V, F):
        """

        :arg V: the mesh's vertices
        :arg F: the mesh's faces

        attributes
        @tri an intrinsic triangulation
        @path a path along triangulation edges
        """
        self.V = V
        self.path = None
        self.tri = F
        # self.intrinsic_edges = self.get_edges(F)

    @staticmethod
    def get_edges(F):
        """
        :arg F: a list of triangular faces, each consisting of 3 vertex indices
        :return: E the corresponding list of edges
        """
        return None

    def set_path(self, path):
        self.path = np.array(path)

    def get_path_polyline(self):
        path_vertices = np.array(
            [self.V[i] for i in self.path]
        )
        path_edges = np.array(
            [[i, i+1] for i in range(len(self.path) - 1)]
        )

        return path_vertices, path_edges

    def get_intrinsic_trinagulation_polyline(self):
        return None, None

    def flip_edge(self, v1, v2):
        pass

