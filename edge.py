import numpy as np
from utils import get_angle, get_side_length


class BaseHalfEdge:
    def __init__(self, twin, origin, vec):
        self.twin = None
        self.set_twin(twin)

        self.origin = origin
        self.next = None

        self.vec = vec
        self.aaaa = "{}->".format(origin)

    def get_origin(self):
        return self.origin

    def get_dst(self):
        return self.next.origin

    def has_endpoint(self, index):
        return self.origin == index or self.get_dst() == index

    def get_vec(self):
        return self.vec

    def set_next(self, next):
        self.next = next
        self.aaaa += str(self.get_dst())

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

    def init_midpoints(self, mesh):
        self.midpoints = []
        if not np.array_equal(self.mesh_face_left, self.mesh_face_right):
            return

        vec = self.vec
        f = self.mesh_face_left
        e = mesh.get_opposite_edge(f, self.origin)

        first_midpoint = e.get_intersection(mesh.V, self.origin, vec)
        self.midpoints.append(first_midpoint)

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

    def print(self, name=""):
        print(name, " ", self.get_origin(), "->", self.get_dst())


class ExtrinsicHalfEdge(BaseHalfEdge):
    def __init__(self, twin, origin, vec):
        super().__init__(twin, origin, vec)

    def get_intersection(self, V, line_start, line_vec):
        # project to plane
        line_vec = np.array(line_vec)
        line_vec /= np.linalg.norm(line_vec)
        line_start = V[line_start]
        plane_normal = line_vec
        plane_point = V[self.origin]

        t = np.dot(plane_normal, plane_point) - np.dot(plane_normal, line_start)
        t /= np.dot(plane_normal, line_vec)
        plane_proj = line_start + line_vec * t

        edge_vec = self.vec
        edge_start = plane_point
        over = plane_proj - edge_start
        s = np.dot(over, edge_vec) / np.dot(edge_vec, edge_vec)
        edge_proj = edge_start + edge_vec * s
        return edge_proj

