from edge import *
import numpy as np
from copy import deepcopy
from utils import turn, is_reflex
import queue

class BaseTriangulation:
    def __init__(self, V):
        self.V = V
        self.in_edges = [[] for _ in range(len(self.V))]
        self.num_f = None

    def add_edge(self, e, dst=None):
        if dst is None:
            dst = e.dst
        self.in_edges[dst].append(e)

    def get_edge(self, origin, dst):
        for e in self.in_edges[dst]:
            if e.origin == origin:
                return e

    def remove_edge_by(self, u, v):
        e = self.get_edge(u, v)
        self.in_edges[v].remove(e)
        self.in_edges[u].remove(e.twin)

    def remove_edge(self, e):
        u, v = e.origin, e.dst
        self.in_edges[v].remove(e)
        self.in_edges[u].remove(e.twin)

    def all_edges(self):
        for v in range(len(self.in_edges)):
            for e in self.in_edges[v]:
                u = e.origin
                if u < v:
                    yield e

    @staticmethod
    def edges_in_face(e, visit_index):
        while not e.was_visited():
            e.mark_visited(visit_index)
            yield e
            e = e.next

    def all_faces(self):
        # mark unvisited
        for edge_list in self.in_edges:
            for e in edge_list:
                e.mark_unvisited()

        face_index = 0
        for edge_list in self.in_edges:
            for e in edge_list:
                if e.was_visited():
                    continue
                yield self.edges_in_face(e, face_index), face_index
                face_index += 1

    def find_shortest_path(self, src, dst):
        num_v = len(self.V)
        d = np.ones(num_v) * np.inf
        visited = np.zeros(num_v, dtype=bool)
        parent = [None for _ in range(num_v)]

        d[src] = 0

        q = queue.PriorityQueue(maxsize=num_v)
        q.put((d[src], src))

        while not q.empty():
            _, v = q.get()
            if visited[v]:
                continue
            visited[v] = True

            for e in self.in_edges[v]:
                u = e.origin
                possible_dist = d[v] + e.length
                if d[u] > possible_dist:
                    d[u] = possible_dist
                    parent[u] = v
                    q.put((d[u], u))

        path = []
        v = dst
        while v != src:
            path.append(v)
            v = parent[v]
            if v is None:
                return None
        path.append(src)
        return self.V[path]


class Triangulation(BaseTriangulation):
    def __init__(self, V, F):
        super().__init__(V)
        self.mesh = ExtrinsicTriangulation(self.V)
        self.insert_mesh_edges(F)
        self.mesh.init_face_angles()

        self.num_f = len(F)
        self.mesh.num_f = self.num_f

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

    def construct_triangle_for_flip(self, twin, prev, next, angle):
        f_left = prev.twin.mesh_face_right
        f_right = prev.next.mesh_face_left

        dir = prev.next.vec

        curr = HalfEdge(twin, prev.dst, f_left, f_right, None)

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

    def flip_by(self, origin, dst):
        e = self.get_edge(origin, dst)
        if e is not None:
            return self.flip(e)

    def flip(self, old_edge):
        triangle_1_prev = old_edge.next
        triangle_2_prev = old_edge.twin.next

        triangle_1_next = triangle_2_prev.next
        triangle_2_next = triangle_1_prev.next

        if triangle_2_next.next != old_edge \
                or triangle_1_next.next != old_edge.twin:
            old_edge.print("Cannot flip", "due to mis-triangulation")
            return None

        triangle_1_angle = old_edge.twin.angle + triangle_1_prev.angle
        triangle_2_angle = old_edge.angle + triangle_2_prev.angle

        if is_reflex(triangle_1_angle) or is_reflex(triangle_2_angle):
            old_edge.print("Cannot flip", "due to reflex angle")
            return None

        self.remove_edge(old_edge)

        e = self.construct_triangle_for_flip(None, triangle_1_prev, triangle_1_next, triangle_1_angle)
        self.construct_triangle_for_flip(e, triangle_2_prev, triangle_2_next, triangle_2_angle)

        e.init_midpoints(self.mesh)

        old_edge.print2("Flipped into", e)
        return e

    def demo_flip(self):
        old = (0, 2)
        new = (1, 3)
        old_edge = self.get_edge(*old)
        self.remove_edge(old_edge)

        # construct triangle 1
        e = self.construct_triangle_for_flip(None, old_edge.next, old_edge.twin.next.next)
        self.construct_triangle_for_flip(e, old_edge.twin.next, old_edge.next.next)

        e.midpoints = []
        e.init_midpoints(self.mesh)
        # e.midpoints.append([-7, 7, 7])

    def get_polyline(self):
        poly_vertices = deepcopy(self.V)
        poly_edges = []
        for e in self.all_edges():
            e_V, e_E = e.get_polyline(index_difference=len(poly_vertices))
            if len(e_V) > 0:
                poly_vertices = np.vstack([poly_vertices, e_V])
            poly_edges.extend(e_E)

        return poly_vertices, np.array(poly_edges)

    def get_faces(self):
        V = deepcopy(self.V)
        F = []
        for f, _ in self.all_faces():
            points = [0]
            for e in f:
                points.append(e.origin)

                midpoints = e.get_midpoints()
                num_midpoints = len(midpoints)
                if num_midpoints > 0:
                    V = np.vstack((V, midpoints))
                    points.extend(list(range(num_midpoints)))

            points[0] = len(points) - 1
            F.append(points)

        return V, np.hstack(F)

    def get_coloring(self, k):
        c = np.zeros(self.num_f, dtype=int) - 1
        for f, f_index in self.all_faces():
            available = np.ones(k, dtype=bool)
            for e in f:
                if not e.twin.was_visited():
                    continue
                neighbour_index = e.twin.get_visit_index()
                neighbour_color = c[neighbour_index]
                available[neighbour_color] = False

            if not np.any(available):
                print("Greedy coloring failed")
                return None

            random_color = np.random.randint(0, k)
            while not available[random_color]:
                random_color = np.random.randint(0, k)

            c[f_index] = random_color

        return c

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
                e.twin = self.get_edge(e.dst, e.origin)

            # add
            self.add_edge(e)

    def init_face_angles(self):
        for e in self.all_edges():
            e.init_face_angle()
