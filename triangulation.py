from edge import *
from exceptions import *
import numpy as np
from copy import deepcopy
from utils import is_reflex_or_flat, get_side_length
from queue import PriorityQueue
from intrinsic_face_triangulator import IntrinsicFaceTriangulator
import view_preferences as prefer
from numeric_error_thresholds import FLAT_ANGLE_THRESHOLD_FOR_EDGE_FLIP

class BaseTriangulation:
    """
    A base class for IntrinsicTriangulation, ExtrinsicTriangulation.

    :field V: the vertices coordinates
    :field in_edges: a linked-list graph of triangulation edges. For each
    vertex v, we maintain a list of its in-edges, i.e. edges whose dst is v.
    """
    def __init__(self, V):
        self.V = V
        self.in_edges = [[] for _ in range(len(self.V))]

    ###################################
    # adding, removing & finding edges
    ###################################
    def add_edge(self, e, dst=None):
        if dst is None:
            dst = e.dst
        self.in_edges[dst].append(e)

    def get_edge(self, origin, dst):
        for e in self.in_edges[dst]:
            if e.origin == origin:
                return e

    def get_all_edges_between(self, origin, dst):
        output = []
        for e in self.in_edges[dst]:
            if e.origin == origin:
                output.append(e)
        return output

    def remove_edge_by(self, origin, dst):
        e = self.get_edge(origin, dst)
        self.in_edges[dst].remove(e)
        self.in_edges[origin].remove(e.twin)

    def remove_edge(self, e):
        u, v = e.origin, e.dst
        self.in_edges[v].remove(e)
        self.in_edges[u].remove(e.twin)

    ############################
    # iterating over all edges
    ############################
    def all_edges(self):
        for v in range(len(self.in_edges)):
            for e in self.in_edges[v]:
                u = e.origin
                if u < v:
                    yield e

    ############################
    # iterating over all faces
    ############################
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

    #################
    #   Dijkstra
    #################
    def find_shortest_path(self, src, dst, parent=None):
        """
        :return: the shortest edge path from src to dst
        If the right edge flips were performed, then this path would be a geodesic;
        Otherwise it is a good approximation for the geodesic.
        :arg parent: (optional) the result of a previously conducted dijkstra search
        """
        if parent is None:
            _, parent = self.dijkstra_distance_and_tree(src)

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

    def dijkstra_distance_and_tree(self, src):
        """
        :arg src: the source vertex from which we begin the Dijkstra search
        :return distance: for each vertex, its distance from src
        :return parent: for each vertex, the previous node on its shortest path from src
        """
        num_v = len(self.V)
        d = np.ones(num_v) * np.inf
        visited = np.zeros(num_v, dtype=bool)
        parent = [None for _ in range(num_v)]

        d[src] = 0

        q = PriorityQueue(maxsize=num_v)
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

        return d, parent


class IntrinsicTriangulation(BaseTriangulation):
    """
    A triangulation over mesh vertices

    :field mesh: the extrinsic triangulation
    :field num_colors: for the coloring of the triangulation, for rendering
    """
    def __init__(self, V, F):
        """
        Input: a mesh, in shared-vertex representation
        Output: the intrinsic triangulation that coincides with the mesh, before any flips
        """
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
                e = IntrinsicHalfEdge(twin, origin, length)
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

    ###############
    #  edge flips
    ###############
    def construct_triangle_for_flip(self, twin, prev, next, angle):
        length = get_side_length(next.length, prev.length, angle)
        curr = IntrinsicHalfEdge(twin, prev.dst, length)

        prev.set_next(curr)
        curr.set_next(next)
        next.set_next(prev)

        prev.corner_angle = angle
        curr.corner_angle = get_angle(prev.length, curr.length, next.length)
        next.corner_angle = get_angle(curr.length, next.length, prev.length)

        such_mesh_edge = self.mesh.get_edge(curr.origin, curr.dst)
        if such_mesh_edge is not None:
            curr.sit_on_mesh_edge(such_mesh_edge)
        else:
            curr.init_near_mesh_edge()

        self.add_edge(curr)

        # stuff for rendering
        self.set_coloring_new_triangle(curr, next, prev)

        return curr

    def flip_by(self, origin, dst):
        e = self.get_edge(origin, dst)

        if e is None:
            raise NonExistentEdgeException(origin, dst)

        return self.flip(e)

    def flip(self, old_edge, flat_angle_threshold=FLAT_ANGLE_THRESHOLD_FOR_EDGE_FLIP):
        """
        edge flip - replacing old_edge by the other diagonal of its quadrilateral

        :arg old_edge: the edge to be flipped
        :arg flat_angle_threshold: to avoid flipping edges if it would create an almost flat angle

        :except ReflexAngleException: cannot flip because the quadrilateral is non-convex
        :except LowDegreeVertexException: cannot flip because the edge has an endpoint of degree 1
        """
        triangle_1_prev = old_edge.next
        triangle_2_prev = old_edge.twin.next

        triangle_1_next = triangle_2_prev.next
        triangle_2_next = triangle_1_prev.next

        if triangle_2_next.next != old_edge \
                or triangle_1_next.next != old_edge.twin:
            raise MistriangulationException(old_edge)

        for v in [old_edge.origin, old_edge.dst]:
            if len(self.in_edges[v]) == 1:
                raise LowDegreeVertexException(old_edge, v, self.in_edges[v])

        triangle_1_angle = old_edge.twin.corner_angle + triangle_1_prev.corner_angle
        triangle_2_angle = old_edge.corner_angle + triangle_2_prev.corner_angle

        if is_reflex_or_flat(triangle_1_angle, flat_angle_threshold) \
                or is_reflex_or_flat(triangle_2_angle, flat_angle_threshold):
            raise ReflexAngleException(old_edge)

        self.remove_edge(old_edge)

        e = self.construct_triangle_for_flip(None, triangle_1_prev, triangle_1_next, triangle_1_angle)
        self.construct_triangle_for_flip(e, triangle_2_prev, triangle_2_next, triangle_2_angle)

        # old_edge.print2("Flipped into", e)
        return e

    def delaunay(self, excluded_edges):
        """
        Turns the triangulation into a Delaunay triangulation, by flipping any illegal edge
        :arg excluded_edges: a set of edges that will not be flipped
        """
        is_any_edge_flipped = False
        edge_generator = self.all_edges()
        more_edges_to_check = []
        while True:
            # get the next edge
            try:
                e = next(edge_generator)
            except StopIteration:
                if len(more_edges_to_check) == 0:
                    break
                edge_generator = (e for e in more_edges_to_check)
                more_edges_to_check = []
                continue

            # skip excluded edges
            if (e.origin, e.dst) in excluded_edges:
                continue
            
            # find the minimal angle in the two adjacent triangles, at the moment
            triangle_1_prev = e.next
            triangle_2_prev = e.twin.next
            triangle_1_next = triangle_2_prev.next
            triangle_2_next = triangle_1_prev.next
            min_actual_angle = min([
                e.corner_angle,
                e.twin.corner_angle,
                triangle_1_prev.corner_angle,
                triangle_2_prev.corner_angle,
                triangle_1_next.corner_angle,
                triangle_2_next.corner_angle,
            ])

            # find the minimal angle in the two adjacent triangles, as it would be after the flip
            triangle_1_angle = e.twin.corner_angle + triangle_1_prev.corner_angle
            triangle_2_angle = e.corner_angle + triangle_2_prev.corner_angle
            new_edge_length = get_side_length(triangle_1_next.length, triangle_1_prev.length, triangle_1_angle)
            new_edge_1_angle = get_angle(triangle_1_prev.length, new_edge_length, triangle_1_next.length)
            next_edge_1_angle = get_angle(new_edge_length, triangle_1_next.length, triangle_1_prev.length)
            new_edge_2_angle = get_angle(triangle_2_prev.length, new_edge_length, triangle_2_next.length)
            next_edge_2_angle = get_angle(new_edge_length, triangle_2_next.length, triangle_2_prev.length)

            min_post_flip_angle = min([
                triangle_1_angle, new_edge_1_angle, next_edge_1_angle,
                triangle_2_angle, new_edge_2_angle, next_edge_2_angle,
            ])

            # flip if illegal
            if min_post_flip_angle > min_actual_angle + 1e-09:
                try:
                    self.flip(e)
                    more_edges_to_check.extend([
                        triangle_1_prev,
                        triangle_1_next,
                        triangle_2_prev,
                        triangle_2_next,
                    ])
                    is_any_edge_flipped = True
                except TriangulationException:
                    pass

        if not is_any_edge_flipped:
            return self.delaunay(excluded_edges)

    ##############
    #  rendering
    ##############
    def get_poly_data(self, mesh, need_extrinsic_faces=True):
        V = deepcopy(self.V)
        E = []
        F = []
        coloring = []

        for f in self.all_faces():
            points = [None]
            is_face_failed = False
            for e in f:
                points.append(e.origin)

                intersections = e.get_intersections(mesh)
                if len(intersections) > 0:
                    index_begin = len(V)
                    V = np.vstack((V, [p.coords for p in intersections]))
                    index_end = len(V)
                    points.extend(list(range(index_begin, index_end)))
                if e.intersections_status == e.Status.FAILED:
                    is_face_failed = True
                    if not prefer.SHOW_FAILED_TRIANGULATION_EDGES:
                        # add the points added so far to E, and start a new one
                        points.extend(points[-2:0:-1])
                        points[0] = len(points) - 1
                        E.append(points)
                        points = [None]

            points[0] = len(points) - 1
            E.append(points)

            if need_extrinsic_faces and not is_face_failed:
                face_color = e.face_color

                if points[0] == 3:
                    F.append(points)
                    coloring.append(face_color)
                else:
                    try:
                        faces, num = self.get_extrinsic_faces(e, intrinsic_face=points[1:])
                        F.append(faces)
                        coloring.extend([face_color] * num)
                    except TriangulationException:
                        pass

        E = np.hstack(E)
        F = np.hstack(F)

        return V, E, F, coloring

    def get_extrinsic_faces(self, e, intrinsic_face):
        return IntrinsicFaceTriangulator(e, intrinsic_face, self.mesh).get_faces()

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

    #############
    # printing
    #############

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
    """
    A triangulation that describes the mesh's shape
    """
    def __init__(self, V):
        super().__init__(V)
        self.avg_face_area = None

    def add_triangle(self, sides):
        """
        Input: an array of 3 IntrinsicHalfEdge objects that coincide to the extrinsic half-edges
        Result: the extrinsic half-edges are added to the extrinsic triangulation
        """
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

        self.avg_face_area = np.mean([next(f).get_triangle_center() for f in self.all_faces()])

    def get_intersection_complete_search(self, line_start, line_vecs, search_source):
        NUM_EDGE_LIMIT = 200
        DISTANCE_LIMIT = 2
        num_edge_count = 0

        trinagle_center = search_source.get_triangle_center()

        for edge_list in self.in_edges:
            for e in edge_list:
                e.mark_unvisited()

        candidates = []
        q = PriorityQueue(maxsize=NUM_EDGE_LIMIT)
        q.put((0, -3, search_source))
        q.put((0, -2, search_source.next))
        q.put((0, -1, search_source.next.next))

        while not q.empty():
            (priority, _, e) = q.get()
            if e.was_visited():
                continue
            twin = e.twin
            e.mark_visited()
            twin.mark_visited()

            for line_vec in line_vecs:
                intersection = e.get_intersection(line_start, line_vec)
                if intersection is not None:
                    intersection.distance_from_face_center = np.linalg.norm(intersection.coords - trinagle_center)
                    candidates.append(intersection)

            next_priority = priority + 1
            if next_priority <= DISTANCE_LIMIT and num_edge_count < NUM_EDGE_LIMIT:
                following_edge = e.next
                while following_edge != twin:
                    q.put((next_priority, num_edge_count, following_edge))
                    num_edge_count += 1
                    following_edge = following_edge.twin.next

        if len(candidates) == 0:
            return None
        print("------------- Successful complete search ------------")
        return min(candidates, key=lambda x: x.distance_from_face_center)



