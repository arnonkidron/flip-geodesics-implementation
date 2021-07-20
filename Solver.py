from FlipEdgeNetwork import FlipEdgeNetwork
from Triangulation import Triangulation
import numpy as np


class Solver:
    def __init__(self, V, F):
        self.net = FlipEdgeNetwork()
        self.tri = Triangulation(F)
        self.V = V
        self.F = F # probably will not be used


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