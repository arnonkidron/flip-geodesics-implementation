from utils import get_angle, get_closest_point, get_angle_between
from math import isclose
from exceptions import *
from edge_intersections import *


class BaseHalfEdge:
    """
    A base class for IntrinsicHalfEdge, ExtrinsicHalfEdge.
    Every edge has twin & next edges as in DCEL. No prev, because all faces are triangles.
    The edge's endpoints are origin and dst.
    :field corner_angle: the angle between the edge and its prev
    :field visited: for DFS search
    :field a_name: just to show the edge's endpoints during debugging
    """
    def __init__(self, twin, origin):
        """
        :param twin: the twin half-edge object, or none if it does not exist
        :param origin: the first endpoint
        """
        self.twin = twin

        self.origin = origin
        self.next = None

        self.corner_angle = None

        self.a_name = "{}->".format(origin)

        self.visited = False
        self.mark_unvisited()

    @property
    def dst(self):
        return self.next.origin

    def set_next(self, next):
        self.next = next
        self.a_name = "{}->{}".format(self.origin, self.dst)

    @property
    def twin(self):
        return self._twin

    @twin.setter
    def twin(self, twin):
        self._twin = twin
        if twin is not None:
            twin._twin = self

    #################
    #   traversal
    #################
    def mark_unvisited(self):
        self.visited = False

    def mark_visited(self):
        self.visited = True

    def was_visited(self):
        return self.visited


class IntrinsicHalfEdge(BaseHalfEdge):
    """
    half-edges for IntrinsicTriangulations
    :field length: the edge's length

    :field near_mesh_edge: the nearest edge from the left in the extrinsic triangulation
    :field angle_from_near_mesh_edge: the angle from it
    Note that this is not how angles are stored in the SignPost datastructure which we did not implement

    :field intersections: the intersection points of this edge with extrinsic edges
    :field intersections_status: whether we finished computing them
    :field face_color: the color to be rendered for the incident intrinsic triangle
    """
    def __init__(self, twin, origin, length):
        """
        :arg length: the edge's length

        __init__ does not initialize next, corner_angle and near_mesh_edge. The user
        must call the appropriate methods, after he finish constructing the other half-edges
        of the triangle in order to complete initialization.
        """
        super().__init__(twin, origin)

        self.length = length

        self.near_mesh_edge = None
        self.angle_from_near_mesh_edge = None

        self.intersections_ = EdgeIntersections(self)

        self.face_color = None

    def to_extrinsic(self, V):
        """
        A constructor for extrinsic edges
        Creates the extrinsic edge that coincides with this intrinsic edge,
        and links this intrinsic edge to it.
        :arg V: the vertices coordinates
        :return: the new extrinsic edge
        """
        ex_half_edge = ExtrinsicHalfEdge(None, self.origin, V[self.origin], V[self.dst])
        ex_half_edge.corner_angle = self.corner_angle

        self.sit_on_mesh_edge(ex_half_edge)

        return ex_half_edge

    ##########################
    # length & corner_angle
    ##########################
    def set_length(self, length):
        if self.length is None:
            self.length = length
            if self.twin is not None:
                self.twin.length = length

    def get_length(self):
        return self.length

    def calc_angles(self):
        next = self.next
        prev = next.next

        prev_len = prev.get_length()
        self_len = self.get_length()
        next_len = next.get_length()

        prev.corner_angle = get_angle(next_len, prev_len, self_len)
        self.corner_angle = get_angle(prev_len, self_len, next_len)
        next.corner_angle = get_angle(self_len, next_len, prev_len)

    #####################
    #   near_mesh_edge
    #####################
    def init_near_mesh_edge(self):
        neighbour = self.next.next.twin
        self.near_mesh_edge = neighbour.near_mesh_edge
        self.angle_from_near_mesh_edge = \
            neighbour.angle_from_near_mesh_edge \
            + self.corner_angle

        next_mesh_edge = self.near_mesh_edge.twin.next
        while self.angle_from_near_mesh_edge > next_mesh_edge.corner_angle:
            self.angle_from_near_mesh_edge -= next_mesh_edge.corner_angle
            self.near_mesh_edge = next_mesh_edge
            next_mesh_edge = self.near_mesh_edge.twin.next

        if isclose(self.angle_from_near_mesh_edge, next_mesh_edge.corner_angle):
            self.sit_on_mesh_edge(next_mesh_edge)

    def sit_on_mesh_edge(self, mesh_edge):
        self.intersections_.set_no_intersection()
        self.near_mesh_edge = mesh_edge
        self.angle_from_near_mesh_edge = 0

    def is_on_mesh_edge(self):
        return self.angle_from_near_mesh_edge == 0

    ###############
    # computing the edge's intersection
    ###############
    def init_intersections(self, mesh):
        return self.intersections_.compute(mesh)

    def init_intersections_one_at_a_time(self, mesh):
        return self.intersections_.compute_one_at_a_time(mesh)

    def get_intersections(self, mesh):
        return self.intersections_.get_intersections(mesh)

    ##############
    #    printing
    ##############
    def print(self, prefix="", suffix=""):
        name = "{}->{}".format(self.origin, self.dst)
        print(prefix, name, suffix)

    def print2(self, middle, other):
        name = "{}->{}".format(self.origin, self.dst)
        other_name = "{}->{}".format(other.origin, other.dst)
        print(name, middle, other_name)

    def get_info(self):
        if self.angle_from_near_mesh_edge != 0:
            n_mid = len(self.get_intersections())
            status = self.intersections_.status.name
            intersection_report = "{} with {} intersections\n".format(status, n_mid)
        else:
            intersection_report = ""

        if self.is_on_mesh_edge():
            mesh_msg = "Angle between faces {:.3f} rad ≈ {:.1f}°\n"\
                .format(self.near_mesh_edge.mesh_face_angle, degrees(self.near_mesh_edge.mesh_face_angle))
        else:
            mesh_msg = ""

        return "Edge {}->{}\n" \
               "Length {:.4f}\n" \
               "Corner Angle {:.3f} rad ≈ {:.1f}°\n" \
               "Near Mesh Edge {}->{}\n" \
               "Angle_from_it {:.3f} rad ≈ {:.1f}°\n" \
            .format(self.origin, self.dst,
                    self.length,
                    self.corner_angle, degrees(self.corner_angle),
                    self.near_mesh_edge.origin, self.near_mesh_edge.dst,
                    self.angle_from_near_mesh_edge, degrees(self.angle_from_near_mesh_edge),
                    ) + intersection_report + mesh_msg


class ExtrinsicHalfEdge(BaseHalfEdge):
    def __init__(self, twin, origin, origin_coords, dst_coords):
        super().__init__(twin, origin)
        self.origin_coords = origin_coords
        self.vec = dst_coords - origin_coords
        self.face_normal = None
        self.triangle_center = None
        self.mesh_face_angle = None

    def get_face_normal(self):
        if self.face_normal is None:
            n = np.cross(self.vec, self.next.vec)
            n /= np.linalg.norm(n)

            self.face_normal = n
            self.next.face_normal = n
            self.next.next.face_normal = n

        return self.face_normal

    def get_triangle_center(self):
        if self.triangle_center is None:
            self.triangle_center = (
                    self.origin_coords
                    + self.next.origin_coords
                    + self.next.next.origin_coords) \
                / 3

        return self.triangle_center

    def get_triangle_area(self):
        return np.linalg.norm(np.cross(self.vec, self.next.vec)) / 2

    def init_face_angle(self):
        if self.twin is None:
            exit("Preprocessing error: Cannot handle open meshes")

        angle = get_angle_between(
            self.get_face_normal(),
            self.twin.get_face_normal()
        )

        self.mesh_face_angle = angle
        self.twin.mesh_face_angle = angle

    def get_intersection(self, line_start, line_vec):
        intersection, distance = get_closest_point(self.origin_coords, self.vec, line_start, line_vec)

        # make sure that it is within the line segment
        diff = intersection - self.origin_coords
        for i in range(3):
            if isclose(self.vec[i], 0):
                continue
            t = diff[i] / self.vec[i]
            if not 0 <= t:
                distance = np.linalg.norm(intersection - self.origin_coords)
            elif not t <= 1:
                distance = np.linalg.norm(intersection - self.next.origin_coords)

        # make sure that it is in front of line_start
        diff = intersection - line_start
        count_dissonances = 0
        for i in range(3):
            if isclose(line_vec[i], 0, abs_tol=1e-02):
                continue
            t = diff[i] / line_vec[i]
            if not(0 <= t):
                count_dissonances += 1
        if count_dissonances >= 2:
            return None

        return Intersection(intersection, self, distance, line_vec)

    def is_point_in_face(self, point, counter=0):
        if counter == 3:
            return False

        if self.origin == point:
            return True

        return self.next.is_point_in_face(point, counter + 1)

