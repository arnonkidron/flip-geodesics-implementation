from edge import *
import numpy as np
from copy import deepcopy
from utils import turn, is_reflex


class BaseTriangulation:
    def __init__(self, V):
        self.V = V
        self.in_edges = [[] for v in range(len(self.V))]

    def add_edge(self, e, dst=None):
        if dst is None:
            dst = e.get_dst()
        self.in_edges[dst].append(e)

    def get_edge(self, origin, dst):
        for e in self.in_edges[dst]:
            if e.origin == origin:
                return e

        return None

    def remove_edge(self, u, v):
        e = self.get_edge(u, v)
        self.in_edges[v].remove(e)
        self.in_edges[u].remove(e.twin)


class Triangulation(BaseTriangulation):
    def __init__(self, V, F):
        super().__init__(V)
        self.mesh = ExtrinsicTriangulation(self.V)
        self.insert_mesh_edges(F)

    def insert_mesh_edges(self, F):
        """
        Initialize the triangulation to include all mesh edges
        Also, add the same edges to self.mesh
        """
        # assume the lists are initialized and empty
        pass

        # fill the lists with the mesh edges
        for f in F:
            sides = []
            for i in range(3):
                origin = f[i]
                dst = f[(i+1) % 3]

                vec = self.V[dst] - self.V[origin]
                twin = self.get_edge(dst, origin)
                e = HalfEdge.construct(twin, origin, f, vec)
                sides.append(e)
                self.add_edge(e, dst)

            # set_next
            for i in range(3):
                sides[i].set_next(sides[(i+1) % 3])

            # initialize angles
            # this must be done only after all edge lengths have been initialized
            sides[0].calc_angles()

            # add to mesh as well
            self.mesh.add_triangle(sides)

    # def get_angle(self, u, v, w):
    #     e1, e2 = None, None
    #     for e in self.in_edges[v]:
    #         if e.has_endpoint(u):
    #             e1 = e
    #         elif e.has_endpoint(w):
    #             e2 = e
    #     if e1 is None or e2 is None:
    #         raise("Cannot compute angle")
    #         # return None
    #
    #     vec1 = e1.get_vec()
    #     vec2 = e2.get_vec()
    #
    #     unit_vec1 = vec1 / np.linalg.norm(vec1)
    #     unit_vec2 = vec2 / np.linalg.norm(vec2)
    #     dot_product = np.dot(unit_vec1, unit_vec2)
    #     angle = np.arccos(dot_product)
    #     return angle

    def construct_triangle_for_flip(self, twin, prev, next, angle):
        f_left = prev.twin.mesh_face_right
        f_right = prev.next.mesh_face_left

        dir = prev.next.vec

        curr = HalfEdge(twin, prev.get_dst(), f_left, f_right, None)

        prev.set_next(curr)
        curr.set_next(next)
        next.set_next(prev)

        prev.angle = angle
        curr.length = get_side_length(next.length, prev.length, prev.angle)
        curr.angle = get_angle(prev.length, curr.length, next.length)
        next.angle = get_angle(curr.length, next.length, prev.length)

        curr.vec = turn(prev.twin.vec, curr.angle, towards=dir)

        self.add_edge(curr)
        return curr

    def flip(self, origin, dst):
        old_edge = self.get_edge(origin, dst)

        triangle_1_prev = old_edge.next
        triangle_1_next = old_edge.twin.next.next
        triangle_1_angle = old_edge.twin.angle + triangle_1_prev.angle

        triangle_2_prev = old_edge.twin.next
        triangle_2_next = old_edge.next.next
        triangle_2_angle = old_edge.angle + triangle_2_prev.angle

        if is_reflex(triangle_1_angle) or is_reflex(triangle_2_angle):
            return None

        self.remove_edge(origin, dst)

        e = self.construct_triangle_for_flip(None, triangle_1_prev, triangle_1_next, triangle_1_angle)
        self.construct_triangle_for_flip(e, triangle_2_prev, triangle_2_next, triangle_2_angle)

        e.midpoints = []
        e.init_midpoints(self.mesh)

        return e

    def demo_flip(self):
        old = (0, 2)
        new = (1, 3)
        old_edge = self.get_edge(*old)
        self.remove_edge(*old)

        # construct triangle 1
        e = self.construct_triangle_for_flip(None, old_edge.next, old_edge.twin.next.next)
        self.construct_triangle_for_flip(e, old_edge.twin.next, old_edge.next.next)

        e.midpoints = []
        e.init_midpoints(self.mesh)
        # e.midpoints.append([-7, 7, 7])

    def all_edges(self):
        for v in range(len(self.in_edges)):
            for e in self.in_edges[v]:
                u = e.get_origin()
                if u < v:
                    yield e

    def get_polyline(self):
        poly_vertices = deepcopy(self.V)
        poly_edges = []
        for e in self.all_edges():
            e_V, e_E = e.get_polyline(index_difference=len(poly_vertices))
            if len(e_V) > 0:
                poly_vertices = np.vstack([poly_vertices, e_V])
            poly_edges.extend(e_E)

        return poly_vertices, np.array(poly_edges)

    def print(self):
        num_V = len(self.in_edges)
        for v in range(num_V):
            for u in range(num_V):
                e = self.get_edge(u, v)
                ch = "V"
                if e is None:
                    ch = "X"
                elif e.twin is None:
                    ch = "A"
                else:
                    ch = "V"
                print(ch, end=" ")
            print("")


class ExtrinsicTriangulation(BaseTriangulation):
    def __init__(self, V):
        super().__init__(V)

    def add_triangle(self, sides):
        edges = [sides[i].to_extrinsic() for i in range(3)]
        for i in range(3):
            e = edges[i]

            # set_next
            e.set_next(edges[(i+1) % 3])

            # set_twin
            if sides[i].twin is not None:
                e.set_twin(self.get_edge(e.get_dst(), e.get_origin()))

            # add
            self.add_edge(e)

    def get_opposite_edge(self, f, v):
        verts = [u for u in f if u != v]
        return self.get_edge(verts[0], verts[1])




