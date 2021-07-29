import numpy as np
from utils import get_angle, get_side_length, rotate, turn, get_closest_point, get_angle_between, is_orientation_counterclockwise, orientation
from math import degrees, isclose, pi
from exceptions import *
from NumericErrorThresholds import *


class BaseHalfEdge:
    def __init__(self, twin, origin):
        self.twin = twin

        self.origin = origin
        self.next = None

        self.corner_angle = None

        self.a_name = "{}->".format(origin)

        self.face_color = None
        self.visited = False
        self.mark_unvisited()

    @property
    def dst(self):
        return self.next.origin

    def set_next(self, next):
        self.next = next
        self.a_name += str(self.dst)

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


class HalfEdge(BaseHalfEdge):
    def __init__(self, twin, origin, length):
        """
        :param twin:
        :param origin:
        :param vec:
        :param mesh_face:

        The user must call set_next, calc_angles after he finish
        constructing the other sides of the triangle.
        """
        super().__init__(twin, origin)

        self.length = length

        self.near_mesh_edge = None
        self.angle_from_near_mesh_edge = None

        # assume the edge lies on a mesh edge, hence it has no intersection points with other mesh edges
        self.num_intersections = None
        self.intersections = []

    def to_extrinsic(self, V):
        ex_half_edge = ExtrinsicHalfEdge(None, self.origin, V[self.origin], V[self.dst])
        ex_half_edge.corner_angle = self.corner_angle

        self.sit_on_mesh_edge(ex_half_edge)

        return ex_half_edge

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
            self.num_intersections = 0
            self.near_mesh_edge = mesh_edge
            self.angle_from_near_mesh_edge = 0

    def is_on_mesh_edge(self):
        return self.angle_from_near_mesh_edge == 0

    ###############
    # computing the edge's intersection
    ###############
    def get_first_segment_vector(self):
        next_mesh_edge = self.near_mesh_edge.twin.next
        return turn(
            self.near_mesh_edge.vec,
            self.angle_from_near_mesh_edge,
            towards=next_mesh_edge.vec
        )

    def get_first_intersecting_mesh_edge(self):
        return self.near_mesh_edge.twin.next.next

    def init_intersections(self, mesh, twin_stuck=False):
        for _ in self.init_intersections_one_at_a_time(mesh, twin_stuck):
            pass

    def init_intersections_one_at_a_time(self, mesh, twin_stuck=False):
        if self.intersections is not None:
            yield None, None; return
        if self.twin.intersections is not None and not twin_stuck:
            yield None, None; return

        self.intersections = []
        dst = self.dst

        if self.near_mesh_edge.is_point_in_face(dst):
            print("--------Superfluous---------")
            self.sit_on_mesh_edge(self.near_mesh_edge)
            yield None, None; return

        # initial values
        prev_intersection = Intersection(self.near_mesh_edge.origin_coords)
        vec = self.get_first_segment_vector()
        e1 = self.get_first_intersecting_mesh_edge()
        e2 = e1

        while len(self.intersections) < 200:
            # try both edges
            candidate1 = e1.get_intersection(prev_intersection.coords, vec)
            candidate2 = e2.get_intersection(prev_intersection.coords, vec)

            if candidate1 is None and candidate2 is None:
                intersection = mesh.get_intersection_complete_search(prev_intersection.coords, vec, e1)
            elif candidate1 is None:
                intersection = candidate2
            elif candidate2 is None:
                intersection = candidate1
            else:
                # if min(candidate1.error, candidate2.error) > INTERSECTION_THRESHOLD:
                #     intersection = mesh.get_intersection_complete_search(prev_intersection.coords, vec, e1)
                if candidate1.error < candidate2.error:
                    intersection = candidate1
                else:
                    intersection = candidate2

            if intersection is None or intersection.error > INTERSECTION_THRESHOLD:
                if twin_stuck:
                    yield None, None; return
                self.twin.init_intersections(mesh, twin_stuck=True)
                other_end = self.twin.intersections
                self.twin.intersections = None
                other_end.reverse()
                self.intersections.extend(other_end)
                yield None, None; return

            self.intersections.append(intersection)
            # yield intersection, self

            e = intersection.mesh_edge.twin
            # if e.is_point_in_face(dst):
            #     yield tmp
            #     break

            prev_intersection = intersection
            vec = rotate(vec, 2 * pi - e.mesh_face_angle, e.vec)
            e1 = e.next
            e2 = e1.next

            intersection.out_vec = vec
            intersection.in_vec = prev_intersection.out_vec

            # n = e.get_face_normal()
            # print(np.dot(vec / np.linalg.norm(vec), n))

            yield intersection, self
            if e.is_point_in_face(dst):
                break

        yield None, None; return

    def get_intersections(self, mesh):
        if self.intersections is None:
            if self.twin.intersections is not None:
                return np.flipud(self.twin.get_intersections(mesh))
            else:
                self.init_intersections(mesh)

        return self.intersections

    def get_polyline(self, index_difference=0):
        intersections = self.get_intersections()

        intersection_indices = list(range(
            index_difference,
            len(intersections) + index_difference
        ))

        verts = [self.origin] \
            + intersection_indices \
            + [self.dst]

        return intersections, [[verts[i], verts[i + 1]] for i in range(len(verts) - 1)]

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

    def print(self, prefix="", suffix=""):
        print(prefix, self.a_name, suffix)

    def print2(self, middle, other):
        print(self.a_name, middle, other.a_name)

    def get_info(self):
        if self.angle_from_near_mesh_edge != 0:
            n_mid = 0 if self.intersections is None else len(self.intersections)
            intersection_report = "Has {} intersections\n".format(n_mid)
        else:
            intersection_report = ""

        if self.is_on_mesh_edge():
            mesh_msg = "Angle between faces {:.3f} rad ≈ {:.1f}°\n"\
                .format(self.near_mesh_edge.mesh_face_angle, degrees(self.near_mesh_edge.mesh_face_angle))
        else:
            mesh_msg = ""

        return "Edge {}->{}\n" \
                "num_intersections {}\n" \
            .format(self.origin, self.dst,
                    # self.length,
                    # self.corner_angle, degrees(self.corner_angle),
                    # self.near_mesh_edge.origin, self.near_mesh_edge.dst,
                    # self.angle_from_near_mesh_edge, degrees(self.angle_from_near_mesh_edge),
                    self.num_intersections
                    ) + intersection_report + mesh_msg


# "Length {:.4f}\n" \
# "Corner Angle {:.3f} rad ≈ {:.1f}°\n" \
# "Near Mesh Edge {}->{}\n" \
# "Angle_from_it {:.3f} rad ≈ {:.1f}°\n" \

class ExtrinsicHalfEdge(BaseHalfEdge):
    def __init__(self, twin, origin, origin_coords, dst_coords):
        super().__init__(twin, origin)
        self.origin_coords = origin_coords
        self.vec = dst_coords - origin_coords
        self.face_normal = None
        self.mesh_face_angle = None

    def get_face_normal(self):
        if self.face_normal is None:
            n = np.cross(self.vec, self.next.vec)
            n /= np.linalg.norm(n)

            self.face_normal = n
            self.next.face_normal = n
            self.next.next.face_normal = n

        return self.face_normal

    def init_face_angle(self):
        angle = get_angle_between(
            self.get_face_normal(),
            self.twin.get_face_normal()
        )

        self.mesh_face_angle = angle
        self.twin.mesh_face_angle = angle

    def get_intersection(self, line_start, line_vec):
        intersection, error = get_closest_point(self.origin_coords, self.vec, line_start, line_vec)

        # make sure that it is within the line segment
        diff = intersection - self.origin_coords
        for i in range(3):
            if isclose(self.vec[i], 0):
                continue
            t = diff[i] / self.vec[i]
            if not 0 <= t:
                error = np.linalg.norm(intersection - self.origin_coords)
            elif not t <= 1:
                error = np.linalg.norm(intersection - self.next.origin_coords)

        # make sure that it is in front of line_start
        diff = intersection - line_start
        count_dissonances = 0
        for i in range(3):
            # if isclose(line_vec[i], 0, abs_tol=1e-02):
            #     continue
            t = diff[i] / line_vec[i]
            if not(0 <= t):
                count_dissonances += 1
        if count_dissonances >= 2:
            return None

        return Intersection(intersection, self, error)

    def is_point_in_face(self, point, counter=0):
        if counter == 3:
            return False

        if self.origin == point:
            return True

        return self.next.is_point_in_face(point, counter + 1)


class Intersection:
    def __init__(self, coords, mesh_edge=None, error=None, out_vec=None):
        self.coords = coords
        self.mesh_edge = mesh_edge
        self.error = error
        self.out_vec = out_vec


