from Triangulation import *
import numpy as np
from utils import is_reflex_or_flat, ROI


def get_shortener(mode, triangulation):
        if mode == ROI.PATH:
            return PathShortener(triangulation)
        elif mode == ROI.LOOP:
            return LoopShortener(triangulation)
        elif mode == ROI.NETWORK:
            return
        elif mode == ROI.SINGLE_SRC:
            return


class BaseShortener:
    def __init__(self, triangulation):
        self.tri = triangulation
        self.is_geodesic = False

        self.length = 0

        self.wedge_angles_forth = []
        self.wedge_angles_back = []

    def flipout(self, a, b, c):
        """
        :arg a, b, c: the wedge that we would like to bypass, to replace with a shorter path a->...->c
        :return bypass: the path that bypasses b
        """
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

        return bypass

    def flipout_the_minimal_wedge_in_path(self):
        pass

    def make_geodesic(self, limit_iterations=None, length_threshold=None):
        if self.is_geodesic:
            return

        initial_length = self.length

        iteration = 0
        while not self.is_geodesic:
            if limit_iterations is not None \
                    and iteration > limit_iterations:
                break
            iteration += 1

            self.flipout_the_minimal_wedge_in_path()

            if length_threshold is not None:
                current_length = self.length
                if current_length / initial_length <= length_threshold:
                    break


class PathShortener(BaseShortener):
    def __init__(self, triangulation):
        super().__init__(triangulation)
        self.path = []
        self.is_geodesic = None

        self.wedge_angles_forth = []
        self.wedge_angles_back = []

    def get_path(self):
        return self.path

    def get_vertex(self, idx):
        return self.path[idx % len(self.path)]

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

        # should be computed lazily?
        self.length = np.sum([self.tri.get_edge(path[i], path[i+1]).length for i in range(len(path) - 1)])

    def update_path(self, b_index, bypass, is_forth):
        # remove b
        b_index = b_index % len(self.path)
        del self.path[b_index]

        # insert bypass to replace it
        if not is_forth:
            bypass.reverse()
            # b_index += 1
        self.path[b_index:b_index] = bypass[1:-1]

        # could be made more efficient
        self.set_path(self.path)

    def flipout_the_minimal_wedge_in_path(self):
        if len(self.path) < 3:
            self.is_geodesic = True
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

        # find the vertices
        b_index = min_index + 1
        if is_forth:
            a, b, c = self.get_vertex(b_index - 1), self.get_vertex(b_index), self.get_vertex(b_index + 1)
        else:
            b_index = -b_index - 1
            c, b, a = self.get_vertex(b_index - 1), self.get_vertex(b_index), self.get_vertex(b_index + 1)

        if is_reflex_or_flat(min_angle):
            self.is_geodesic = True
            return self.path

        bypass = self.flipout(a, b, c)
        self.update_path(b_index, bypass, is_forth)

        return self.path


    ################
    # move to SingleSourceShortener
    ################
    def set_tree(self, src, parent):
        self.path = []
        self.is_geodesic = False

        num_v = len(self.tri.V)

        self.wedge_angles_forth = [
            (v,
                self.tri.get_wedge_angle(
                    v,
                    parent[v],
                    parent[parent[v]],
                )
            ) for v in range(num_v) if src not in [v, parent[v]]
        ]
        self.wedge_angles_back = [
            (v,
                self.tri.get_wedge_angle(
                    parent[parent[v]],
                    parent[v],
                    v,
                )
            ) for v in range(num_v) if src not in [v, parent[v]]
        ]

    def update_tree(self, bypass, is_forth):
        if not is_forth:
            bypass.reverse()

        pass

    def flipout_the_minimal_wedge_in_tree(self, parent):
        if self.is_geodesic:
            return

        # find minimal wedge
        min_vertex_forth, min_angle_forth = min(self.wedge_angles_forth)
        min_vertex_back, min_angle_back = min(self.wedge_angles_back)

        if min_angle_forth < min_angle_back:
            min_vertex, min_angle, is_forth = min_vertex_forth, min_angle_forth, True
        else:
            min_vertex, min_angle, is_forth = min_vertex_back, min_angle_back, False

        if is_reflex_or_flat(min_angle):
            self.is_geodesic = True
            return self.path

        # find the vertices
        if is_forth:
            a, b, c = min_vertex, parent[min_vertex], parent[parent[min_vertex]]
        else:
            c, b, a = min_vertex, parent[min_vertex], parent[parent[min_vertex]]

        bypass = self.flipout(a, b, c)

        # update tree
        self.update_tree(bypass, is_forth)

        return self.path

    def make_single_source_geodesic(self, src, limit_iterations=None, length_threshold=None):
        _, parent = self.tri.dijkstra_distance_and_tree(src)
        self.set_tree(src, parent)

        iteration = 0
        while not self.is_geodesic:
            if limit_iterations is not None \
                    and iteration > limit_iterations:
                break
            iteration += 1

            self.flipout_the_minimal_wedge_in_tree(parent)

            _, parent = self.tri.dijkstra_distance_and_tree(src)
            self.set_tree(src, parent)

        print("Done {} iterations".format(limit_iterations))


class LoopShortener(PathShortener):
    def __init__(self, triangulation):
        super().__init__(triangulation)

    def get_path(self):
        if self.path[0] != self.path[-1]:
            self.path.append(self.path[0])
        return self.path

    def set_path(self, path):
        if path[-1] == path[0]:
            path.pop()
        self.path = path
        self.is_geodesic = len(path) <= 1

        self.wedge_angles_forth = [
            self.tri.get_wedge_angle(
                self.get_vertex(i-1),
                self.get_vertex(i),
                self.get_vertex(i+1),
            ) for i in range(1, len(path) - 1 + 2)
        ]
        self.wedge_angles_back = [
            self.tri.get_wedge_angle(
                self.get_vertex(i+1),
                self.get_vertex(i),
                self.get_vertex(i-1),
            ) for i in range(len(path) - 2, 0 - 2, -1)
        ]
        self.length = np.sum([self.tri.get_edge(path[i], path[i+1]).length for i in range(-1, len(path) - 1)])

    def flipout_the_minimal_wedge_in_path(self):
        super().flipout_the_minimal_wedge_in_path()
        if len(self.path) == 2 and len(self.tri.get_all_edges_between(self.path[0], self.path[1])) == 1:
            self.path.pop()

