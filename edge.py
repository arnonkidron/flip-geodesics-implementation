import numpy as np
from utils import get_angle, get_side_length, rotate, turn, get_closest_point, get_angle_between, is_orientation_counterclockwise, orientation
from math import degrees, isclose


class BaseHalfEdge:
    def __init__(self, twin, origin, vec):
        self.twin = twin

        self.origin = origin
        self.next = None

        self.angle = None

        self.vec = None
        if vec is not None:
            self.set_vec(vec)

        self.a_name = "{}->".format(origin)

        self.face_color = None
        self.visited = False
        self.mark_unvisited()

    @property
    def dst(self):
        return self.next.origin

    def set_vec(self, vec):
        self.vec = vec / np.linalg.norm(vec)

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
    def __init__(self, twin, origin, mesh_face_left, mesh_face_right, vec):
        """
        :param twin:
        :param origin:
        :param vec:
        :param mesh_face:

        The user must call set_next, calc_angles after he finish
        constructing the other sides of the triangle.
        """
        super().__init__(twin, origin, vec)

        self.mesh_face_left = mesh_face_left
        self.mesh_face_right = mesh_face_right
        self.length = None
        self.midpoints = None

        self.near_mesh_edge = None
        self.angle_from_near_mesh_edge = None

        self.is_done_init_midpoints = True
        self.midpoints_generator = None

    @staticmethod
    def construct(twin, origin, mesh_face_left, vec):
        """
        a constructor for half-edges that lie on mesh edges
        """
        obj = HalfEdge(twin, origin, mesh_face_left, None, vec)

        if twin is not None:
            obj.mesh_face_right = twin.mesh_face_left
            twin.mesh_face_right = obj.mesh_face_left

        obj.set_length(np.linalg.norm(vec))
        obj.midpoints = []

        return obj

    def to_extrinsic(self):
        ex_half_edge = ExtrinsicHalfEdge(None, self.origin, self.vec)
        ex_half_edge.angle = self.angle
        self.near_mesh_edge = ex_half_edge
        self.angle_from_near_mesh_edge = 0
        return ex_half_edge

    def init_near_mesh_edge(self):
        neighbour = self.next.next.twin
        self.near_mesh_edge = neighbour.near_mesh_edge
        self.angle_from_near_mesh_edge = \
            neighbour.angle_from_near_mesh_edge \
            + self.angle

        next_mesh_edge = self.near_mesh_edge.twin.next
        while self.angle_from_near_mesh_edge > next_mesh_edge.angle:
            self.angle_from_near_mesh_edge -= next_mesh_edge.angle
            self.near_mesh_edge = next_mesh_edge
            next_mesh_edge = self.near_mesh_edge.twin.next

        if isclose(self.angle_from_near_mesh_edge, next_mesh_edge.angle):
            self.sit_on_mesh_edge(next_mesh_edge)

    def sit_on_mesh_edge(self, mesh_edge):
            self.near_mesh_edge = mesh_edge
            self.angle_from_near_mesh_edge = 0

    def get_vec_by_near_mesh_edge(self):
        next_mesh_edge = self.near_mesh_edge.twin.next
        return turn(
            self.near_mesh_edge.vec,
            self.angle_from_near_mesh_edge,
            towards=next_mesh_edge.vec
        )

    # auxiliary methods for init_midpoints
    def get_initial_target_edge(self, mesh):
        return self.near_mesh_edge.twin.next.next
        # f = self.mesh_face_left
        # i = np.nonzero(f == self.origin)[0][0]
        # return mesh.get_edge(f[(i+1) % 3], f[(i+2) % 3])

    @staticmethod
    def get_next_target_edge(prev_midpoint, vec, e, mesh):
        e1 = e.next
        e2 = e1.next
        opposite_vertex = e2.origin
        n = e.get_face_normal()

        opposite_vertex = mesh.V[opposite_vertex]
        e1_origin = mesh.V[e1.origin]

        orientation_with_e1 = is_orientation_counterclockwise(e1_origin, prev_midpoint, opposite_vertex, n)
        orientation_with_vec = is_orientation_counterclockwise(e1_origin, prev_midpoint, prev_midpoint + vec, n)

        if orientation_with_e1 != orientation_with_vec:
            return e1
        else:
            return e2


        n = e.get_face_normal()
        orientation_e = is_orientation_counterclockwise(opposite_vertex - e1.vec, opposite_vertex, opposite_vertex + e2.vec, n)
        orientation_e2 = is_orientation_counterclockwise(opposite_vertex + e2.vec, opposite_vertex, opposite_vertex - e1.vec, n)
        if orientation_e == orientation_e2:
            print("totally wrong!")

        orientation_v = is_orientation_counterclockwise(opposite_vertex, prev_midpoint, prev_midpoint + vec, n)

        if orientation_e == orientation_v:
            return e1
        else:
            return e2

    def restore_mesh_edge(self, mesh):
        """
        A method to be called when mesh edges are created by edge flips.
        We replace its data with the exact data from the mesh.
        """
        print("Shouldn't have gotten here!")

        f = self.next.mesh_face_left
        self.mesh_face_left = f
        self.twin.mesh_face_right = f

        f = self.twin.next.mesh_face_left
        self.mesh_face_right = f
        self.twin.mesh_face_left = f

        self.vec = mesh.get_edge(self.origin, self.dst).vec
        self.twin.vec = -self.vec

    def init_midpoints(self, mesh):
        self.is_done_init_midpoints = True
        if self.midpoints is not None \
                or self.twin.midpoints is not None:
            return

        self.midpoints = []
        dst = self.dst
        if dst in self.mesh_face_left:
            print("Oh!")
            self.sit_on_mesh_edge("If anyone gets here, that's fine, but we'd need to code that")
            return self.near_mesh_edge

        self.is_done_init_midpoints = False

        # initial values
        prev_midpoint = mesh.V[self.origin]
        vec = self.vec
        e = self.get_initial_target_edge(mesh)

        while len(self.midpoints) < 200:
            midpoint = e.get_intersection(mesh.V, prev_midpoint, vec)
            if midpoint is None:
                e = e.next
                midpoint = e.get_intersection(mesh.V, prev_midpoint, vec)
                if midpoint is None:
                    print("Midpoint calculation has failed")
                    return
            self.midpoints.append(midpoint)
            yield e, midpoint

            e = e.twin
            if e.is_point_in_face(dst):
                break

            prev_midpoint = midpoint
            vec = rotate(vec, e.mesh_face_angle, e.vec)
            e = e.next

        self.is_done_init_midpoints = True
        yield e, midpoint

    def get_midpoints(self):
        if self.midpoints is None and self.twin.midpoints is not None:
            return np.flipud(self.twin.get_midpoints())

        return self.midpoints

    def get_polyline(self, index_difference=0):
        midpoints = self.get_midpoints()

        midpoint_indices = list(range(
            index_difference,
            len(midpoints) + index_difference
        ))

        verts = [self.origin] \
            + midpoint_indices \
            + [self.dst]

        return midpoints, [[verts[i], verts[i + 1]] for i in range(len(verts) - 1)]

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

        prev.angle = get_angle(next_len, prev_len, self_len)
        self.angle = get_angle(prev_len, self_len, next_len)
        next.angle = get_angle(self_len, next_len, prev_len)

    def print(self, prefix="", suffix=""):
        print(prefix, self.a_name, suffix)

    def print2(self, middle, other):
        print(self.a_name, middle, other.a_name)

    def get_info(self):

        n_mid = 0 if self.midpoints is None else len(self.midpoints)
        if self.is_done_init_midpoints:
            midpoint_report = "Done with {} midpoints\n".format(n_mid)
        else:
            midpoint_report = "Ongoing with {} midpoints\n".format(n_mid)

        if self.angle_from_near_mesh_edge == 0:
            mesh_msg = "Angle between faces {:.3f} rad ≈ {:.1f}°\n"\
                .format(self.near_mesh_edge.mesh_face_angle, degrees(self.near_mesh_edge.mesh_face_angle))
        else:
            mesh_msg = ""

        return "Edge {}->{}\n" \
               "Length {:.4f}\n" \
               "Left_Angle {:.3f} rad ≈ {:.1f}°\n" \
               "Start_Vec ({:.4f},{:.4f},{:.4f})\n" \
               "Near Mesh Edge {}->{}\n" \
               "Angle_from_it {:.3f} rad ≈ {:.1f}°\n" \
               "Mesh_face_left: [{},{},{}]\n" \
               "Mesh_face_right: [{},{},{}]\n" \
            .format(self.origin, self.dst,
                    self.length,
                    self.angle, degrees(self.angle),
                    self.vec[0], self.vec[1], self.vec[2],
                    self.near_mesh_edge.origin, self.near_mesh_edge.dst,
                    self.angle_from_near_mesh_edge, degrees(self.angle_from_near_mesh_edge),
                    self.mesh_face_left[0], self.mesh_face_left[1], self.mesh_face_left[2],
                    self.mesh_face_right[0], self.mesh_face_right[1], self.mesh_face_right[2]
                    ) + midpoint_report + mesh_msg


class ExtrinsicHalfEdge(BaseHalfEdge):
    def __init__(self, twin, origin, vec):
        super().__init__(twin, origin, vec)
        self.face_normal = None
        self.mesh_face_angle = None

    def get_face_normal(self):
        if self.face_normal is None:
            # we may assume all vecs are normalized, because set_vec takes care of that
            n = np.cross(self.vec, self.next.vec)
            n /= np.linalg.norm(n)

            self.face_normal = n
            self.next.face_normal = n
            self.next.next.face_normal = n

        return self.face_normal

    def init_face_angle(self):
        # TODO: find why this fails so terribly
        angle = get_angle_between(
            self.get_face_normal(),
            self.twin.get_face_normal()
        )

        self.mesh_face_angle = angle
        self.twin.mesh_face_angle = angle

    def get_intersection(self, V, line_start, line_vec):
        self_start = V[self.origin]
        intersection = get_closest_point(self_start, self.vec, line_start, line_vec)

        # make sure that it is within the line segment
        self_end = V[self.dst]
        unnormalized_vec = self_end - self_start
        diff = intersection - self_start
        for i in range(3):
            if isclose(unnormalized_vec[i], 0):
                continue
            t = diff[i] / unnormalized_vec[i]
            if not(0 <= t <= 1):
                return None

        return intersection

    def is_point_in_face(self, point, counter=0):
        if counter == 3:
            return False

        if self.origin == point:
            return True

        return self.next.is_point_in_face(point, counter + 1)

