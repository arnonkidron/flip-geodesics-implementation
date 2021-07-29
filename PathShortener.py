from exceptions import *
from Triangulation import *
import numpy as np
from utils import is_reflex_or_flat


class PathShortener:
    def __init__(self, triangulation):
        self.tri = triangulation
        self.path = []
        self.is_geodesic = None

        self.length = 0
        self.wedge_angles_forth = []
        self.wedge_angles_back = []

    def get_path(self):
        return self.path

    def set_path(self, path):
        self.path = path
        self.is_geodesic = len(path) <= 2

        self.wedge_angles_forth = [
            self.tri.get_wedge_angle(
                self.path[i-1],
                self.path[i],
                self.path[i+1],
            ) for i in range(1, len(path) - 1)
        ]
        self.wedge_angles_back = [
            self.tri.get_wedge_angle(
                self.path[i+1],
                self.path[i],
                self.path[i-1],
            ) for i in range(len(path) - 2, 0, -1)
        ]
        self.length = np.sum([self.tri.get_edge(path[i], path[i+1]).length for i in range(len(path) - 1)])

    def update_path(self, b_index, is_forth, bypass):
        # remove b
        del self.path[b_index]

        # insert bypass to replace it
        if not is_forth:
            bypass.reverse()
            b_index += 1
        self.path[b_index:b_index] = bypass[1:-1]

        self.set_path(self.path)

    def flipout(self, b_index, is_forth):
        if is_forth:
            a, b, c = self.path[b_index - 1], self.path[b_index], self.path[b_index + 1]
        else:
            b_index = -b_index - 1
            c, b, a = self.path[b_index - 1], self.path[b_index], self.path[b_index + 1]

        wedge_angle = self.tri.get_wedge_angle(a, b, c)
        if is_reflex_or_flat(wedge_angle):
            raise WedgeReflexAngleException(wedge_angle)

        # find edges of this wedge
        e1 = self.tri.get_edge(a, b)
        e2 = self.tri.get_edge(b, c)
        if e1 is None or e2 is None:
            return

        if e1 == e2 or e1 == e2.twin:
            raise SelfEdgeException(a, b, c)

        # calculate bypassing path
        bypass = [a]
        e = e1.next
        while e.dst != c:
            bypass.append(e.dst)
            e = e.twin.next
        bypass.append(c)

        i = 1
        while i != len(bypass) - 1:
            to_flip_edge = self.tri.get_edge(bypass[i], b)

            try:
                e = self.tri.flip(to_flip_edge)
            except ReflexAngleException:
                i = i + 1
                continue
            except TriangulationException as err:
                raise err

            del bypass[i]
            if i > 1:
                i = i - 1

        # update self.path
        self.update_path(b_index, is_forth, bypass)

    def flipout_the_minimal_wedge(self):
        if len(self.path) < 3:
            return self.path

        # find minimal wedge
        min_index_forth = np.argmin(self.wedge_angles_forth)
        min_index_back = np.argmin(self.wedge_angles_back)
        
        min_angle_forth = self.wedge_angles_forth[min_index_forth]
        min_angle_back = self.wedge_angles_back[min_index_back]

        if min_angle_forth < min_angle_back:
            min_index, min_angle, is_forth = min_index_forth, min_angle_forth, True
        else:
            min_index, min_angle, is_forth = min_index_back, min_angle_back, False

        # b_index is the index of the middle vertex of the wedge
        b_index = min_index + 1

        if is_reflex_or_flat(min_angle):
            self.is_geodesic = True
            return self.path

        self.flipout(b_index, is_forth)

        return self.path

    def make_geodesic(self, limit_iterations=None, length_threshold=None):
        initial_length = self.length

        iteration = 0
        while not self.is_geodesic:
            if limit_iterations is not None \
                    and iteration < limit_iterations:
                break
            iteration += 1

            self.flipout_the_minimal_wedge()

            if length_threshold is not None:
                current_length = self.length
                if current_length / initial_length <= length_threshold:
                    break

        return self.path

    def is_geodesic(self):
        return False
