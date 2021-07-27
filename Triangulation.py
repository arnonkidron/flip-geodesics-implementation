from edge import *
from exceptions import *
import numpy as np
from copy import deepcopy
from utils import turn, is_reflex
import queue


class BaseTriangulation:
    def __init__(self, V):
        self.V = V
        self.in_edges = [[] for _ in range(len(self.V))]

    def get_wedge_angle(self, a, b, c):
        """
        from left, i.e. the clockwise angle
        """
        e1 = self.get_edge(a, b)
        e2 = self.get_edge(b, c)

        if e1 is None or e2 is None:
            print("The vertices {}->{}->{} do not form a 2-path".format(a, b, c))
            return None

        if e1 == e2:
            return 0

        e = e1.next
        sum = 0
        while e != e2:
            sum += e.corner_angle
            e = e.twin.next
        sum += e.corner_angle

        return sum

    def get_wedge_angle_counterclockwise(self, a, b, c):
        return self.get_wedge_angle(c, b, a)

    def add_edge(self, e, dst=None):
        if dst is None:
            dst = e.dst
        self.in_edges[dst].append(e)

    def get_edge(self, origin, dst):
        for e in self.in_edges[dst]:
            if e.origin == origin:
                return e

    def remove_edge_by(self, origin, dst):
        e = self.get_edge(origin, dst)
        self.in_edges[dst].remove(e)
        self.in_edges[origin].remove(e.twin)

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
    def edges_in_face(e):
        while not e.was_visited():
            e.mark_visited()
            yield e
            e = e.next

    def all_faces(self):
        # mark all unvisited
        for edge_list in self.in_edges:
            for e in edge_list:
                e.mark_unvisited()

        # generate
        for edge_list in self.in_edges:
            for e in edge_list:
                if e.was_visited():
                    continue
                yield self.edges_in_face(e)

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
        path.reverse()
        return path


class Triangulation(BaseTriangulation):
    def __init__(self, V, F):
        super().__init__(V)
        self.mesh = ExtrinsicTriangulation(self.V)
        self.insert_mesh_edges(F)
        self.mesh.init_face_angles()

        self.num_colors = None

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

                length = np.linalg.norm(self.V[dst] - self.V[origin])
                twin = self.get_edge(dst, origin)
                e = HalfEdge(twin, origin, length)
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
        length = get_side_length(next.length, prev.length, angle)
        curr = HalfEdge(twin, prev.dst, length)

        prev.set_next(curr)
        curr.set_next(next)
        next.set_next(prev)

        prev.corner_angle = angle
        curr.corner_angle = get_angle(prev.length, curr.length, next.length)
        next.corner_angle = get_angle(curr.length, next.length, prev.length)

        curr.init_near_mesh_edge()

        self.add_edge(curr)

        # stuff for rendering
        curr.intersections = None
        self.set_coloring_new_triangle(curr, next, prev)

        return curr

    def flip_by(self, origin, dst):
        e = self.get_edge(origin, dst)

        if e is None:
            raise NonExistentEdgeException(e)

        return self.flip(e)

    def flip(self, old_edge):
        triangle_1_prev = old_edge.next
        triangle_2_prev = old_edge.twin.next

        triangle_1_next = triangle_2_prev.next
        triangle_2_next = triangle_1_prev.next

        if triangle_2_next.next != old_edge \
                or triangle_1_next.next != old_edge.twin:
            raise MistriangulationException(old_edge)

        triangle_1_angle = old_edge.twin.corner_angle + triangle_1_prev.corner_angle
        triangle_2_angle = old_edge.corner_angle + triangle_2_prev.corner_angle

        if is_reflex(triangle_1_angle) or is_reflex(triangle_2_angle):
            raise ReflexAngleException(old_edge)

        self.remove_edge(old_edge)

        e = self.construct_triangle_for_flip(None, triangle_1_prev, triangle_1_next, triangle_1_angle)
        self.construct_triangle_for_flip(e, triangle_2_prev, triangle_2_next, triangle_2_angle)

        old_edge.print2("Flipped into", e)
        return e

    def get_polyline(self):
        poly_vertices = deepcopy(self.V)
        poly_edges = []
        for e in self.all_edges():
            e_V, e_E = e.get_polyline(index_difference=len(poly_vertices))
            if len(e_V) > 0:
                poly_vertices = np.vstack([poly_vertices, e_V])
            poly_edges.extend(e_E)

        return poly_vertices, np.array(poly_edges)

    def get_poly_data(self):
        V = deepcopy(self.V)
        E = []
        coloring = []
        for f in self.all_faces():
            points = [None]
            for e in f:
                points.append(e.origin)

                intersections = e.get_intersections()
                if intersections is not None and len(intersections) > 0:
                    index_begin = len(V)
                    V = np.vstack((V, intersections))
                    index_end = len(V)
                    points.extend(list(range(index_begin, index_end)))

            points[0] = len(points) - 1
            E.append(points)

            face_color = e.face_color
            coloring.append(face_color)



        E = np.hstack(E)

        return V, E, E, coloring

    def init_coloring(self, num_faces, num_colors):
        self.face_coloring = np.zeros(num_faces, dtype=int) - 1
        self.num_colors = num_colors
        for face in self.all_faces():
            self.set_coloring(face)

    def set_coloring(self, face):
        available = np.ones(self.num_colors, dtype=bool)
        for e in face:
            if not e.twin.was_visited():
                continue
            neighbour_color = e.twin.face_color
            available[neighbour_color] = False

        if not np.any(available):
            print("Greedy coloring failed. Need more colors")
            return None

        available, = np.where(available)
        random_color = available[np.random.randint(0, len(available))]

        e.face_color = random_color
        e.next.face_color = random_color
        e.next.next.face_color = random_color

    def set_coloring_new_triangle(self, new_edge, next, nextnext):
        if next.face_color is not None:
            next.mark_unvisited()
            nextnext.mark_unvisited()
            if new_edge.twin is None:
                new_edge.mark_visited()
            self.set_coloring(self.edges_in_face(next))

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

    def get_vertex_info(self, idx):
        coords = self.V[idx]
        deg = len(self.in_edges[idx])
        return "Vertex #{}\n" \
               "({:.4f},{:.4f},{:.4f})\n" \
               "Degree {}\n"\
            .format(idx, coords[0], coords[1], coords[2], deg)


class ExtrinsicTriangulation(BaseTriangulation):
    def __init__(self, V):
        super().__init__(V)

    def add_triangle(self, sides):
        edges = [sides[i].to_extrinsic(self.V) for i in range(3)]
        for i in range(3):
            e = edges[i]

            # set_next
            e.set_next(edges[(i+1) % 3])

            # set twin
            if sides[i].twin is not None:
                e.twin = self.get_edge(e.dst, e.origin)

            # add
            self.add_edge(e)

    def init_face_angles(self):
        for e in self.all_edges():
            e.init_face_angle()
