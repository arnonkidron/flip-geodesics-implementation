import numpy as np
from utils import get_angle, get_side_length, rotate, get_closest_point, get_angle_between


class BaseHalfEdge:
    def __init__(self, twin, origin, vec):
        self.twin = None
        self.set_twin(twin)

        self.origin = origin
        self.next = None

        self.vec = None
        if vec is not None:
            self.set_vec(vec)

        self.a_name = "{}->".format(origin)

    def get_origin(self):
        return self.origin

    def get_dst(self):
        return self.next.origin

    def has_endpoint(self, index):
        return self.origin == index or self.get_dst() == index

    def get_vec(self):
        return self.vec

    def set_vec(self, vec):
        self.vec = vec / np.linalg.norm(vec)

    def set_next(self, next):
        self.next = next
        self.a_name += str(self.get_dst())

    def set_twin(self, twin):
        self.twin = twin
        if twin is not None:
            twin.twin = self


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
        self.angle = None
        self.length = None
        self.midpoints = None

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
        return ExtrinsicHalfEdge(None, self.origin, self.vec)

    def get_angle(self):
        return self.angle

    # auxiliary methods for init_midpoints
    def get_initial_target_edge(self, mesh):
        return mesh.get_opposite_edge(self.mesh_face_left, self.origin)

    def get_next_target_edge(self, prev_midpoint, vec, e):
        e1 = e.next
        e2 = e1.next
        opposite_vertex = e2.origin

        segment = opposite_vertex - prev_midpoint
        segment /= np.linalg.norm(segment)
        angle = get_angle_between(segment, vec)
        if angle < 0:
            return e1
        else:
            return e2



    def init_midpoints(self, mesh):
        self.midpoints = []
        if not np.array_equal(self.mesh_face_left, self.mesh_face_right):
            return

        dst = self.get_dst()

        # initial values
        prev_midpoint = mesh.V[self.origin]
        vec = self.vec
        e = self.get_initial_target_edge(mesh)

        while True:
            midpoint = e.get_intersection(mesh.V, prev_midpoint, vec)
            self.midpoints.append(midpoint)

            if e.twin.is_point_in_face(dst):
                break

            prev_midpoint = midpoint
            vec = rotate(vec, e.mesh_face_angle, e.vec)
            e = self.get_next_target_edge(prev_midpoint, vec, e)

        return

        vec = self.vec
        f = self.mesh_face_left
        e = mesh.get_opposite_edge(f, self.origin)

        first_midpoint = e.get_intersection(mesh.V, mesh.V[self.origin], vec)
        self.midpoints.append(first_midpoint)

        dst = self.get_dst()

        if e.twin.is_point_in_face(dst):
            return

        vec = rotate(vec, e.mesh_face_angle, e.vec)

        e1 = e.next
        e2 = e1.next

        e = e2   # by oracle
        second_midpoint = e.get_intersection(mesh.V, first_midpoint, vec)
        self.midpoints.append(second_midpoint)

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
            + [self.get_dst()]

        return midpoints, [[verts[i], verts[i + 1]] for i in range(len(verts) - 1)]

    def set_length(self, length):
        if self.length is None:
            self.length = length
            if self.twin is not None:
                self.twin.length = length

    def get_length(self):
        if self.length is None:
            # uninitialized.
            # we may assume that the two other sides are initialized,
            # because edge flips add only one new edge.
            a = self.next
            b = a.next
            self.length = get_side_length(a.length, b.length, b.angle)

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

    def print(self, prefix=""):
        print(prefix, self.a_name)

    def print2(self, middle, other):
        print(self.a_name, middle, other.a_name)


class ExtrinsicHalfEdge(BaseHalfEdge):
    def __init__(self, twin, origin, vec):
        super().__init__(twin, origin, vec)
        self.face_normal = None
        self.mesh_face_angle = None

    def get_face_normal(self):
        if self.face_normal is None:
            # we may assume all vecs are normalized, because set_vec takes care of that
            n = np.cross(self.vec, self.next.vec)

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

    def get_intersection(self, V, line_start, line_vec):
        self_start = V[self.origin]
        return get_closest_point(self_start, self.vec, line_start, line_vec)

    def is_point_in_face(self, point, counter=0):
        if counter == 3:
            return False

        if self.origin == point:
            return True

        return self.next.is_point_in_face(point, counter + 1)

