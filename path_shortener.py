from triangulation import *
import numpy as np
from utils import is_reflex_or_flat, ROI
from numeric_error_thresholds import *
from joint import Joint


def get_shortener(mode, triangulation):
        if mode == ROI.PATH:
            return PathShortener(triangulation)
        elif mode == ROI.LOOP:
            return LoopShortener(triangulation)
        elif mode == ROI.NETWORK:
            return MultiplePathShortener(triangulation)


class BaseShortener:
    def __init__(self, triangulation):
        self.tri = triangulation
        self.is_geodesic = False

        self.joints = []

        self.reflex_angle_threshold_for_edge_flip = FLAT_ANGLE_THRESHOLD_FOR_EDGE_FLIP
        self.reflex_angle_threshold_for_flip_out = FLAT_ANGLE_THRESHOLD_FOR_FLIP_OUT

        self.keep_list_of_flipped_edges = False
        self.flipped_edges = []

    @property
    def length(self):
        return 0

    def flipout(self, joint):
        """
        :arg joint: the joint a->b->c that we would like to bypass, to replace with a shorter path a->...->c
        :return bypass: the path that bypasses b
                the triangulation is modified, so it includes the bypass and it no longer includes the joint
        """
        a, b, c = joint.vertices

        if is_reflex_or_flat(joint.wedge_angle, self.reflex_angle_threshold_for_flip_out):
            raise WedgeReflexAngleException(joint.wedge_angle)

        if joint.e1 == joint.e2 or joint.e1 == joint.e2.twin:
            raise SelfEdgeException(a, b, c)

        # calculate bypassing path, by going along the wedge edges
        bypass = [a]
        wedge_edges = [None]

        e = joint.e1.next
        while e.dst != c:
            bypass.append(e.dst)
            wedge_edges.append(e)
            e = e.twin.next
        bypass.append(c)

        i = 1
        while i != len(bypass) - 1:
            to_flip_edge = wedge_edges[i]

            try:
                self.tri.flip(to_flip_edge, flat_angle_threshold=self.reflex_angle_threshold_for_edge_flip)
                if self.keep_list_of_flipped_edges:
                    self.flipped_edges.append(to_flip_edge)
            except ReflexAngleException:
                i = i + 1
                continue
            except TriangulationException as err:
                raise err

            del bypass[i]
            del wedge_edges[i]
            if i > 1:
                i = i - 1

        return bypass

    def flipout_the_minimal_wedge(self):
        pass

    def make_geodesic(self, limit_iterations=None, length_threshold=None):
        if self.is_geodesic:
            return

        initial_length = self.length

        self.flipped_edges = []
        iteration = 0
        while not self.is_geodesic:
            if limit_iterations is not None \
                    and iteration > limit_iterations:
                break
            iteration += 1

            self.flipout_the_minimal_wedge()

            if length_threshold is not None:
                current_length = self.length
                if current_length / initial_length <= length_threshold:
                    break


class PathShortener(BaseShortener):
    def __init__(self, triangulation):
        super().__init__(triangulation)
        self.path = []
        self.is_geodesic = None

    def get_path(self):
        return self.path

    def get_vertex(self, idx):
        return self.path[idx % len(self.path)]

    @property
    def length(self):
        return np.sum([self.tri.get_edge(self.path[i], self.path[i+1]).length for i in range(len(self.path) - 1)])

    def is_path_too_short(self):
        return len(self.path) <= 2

    def set_path(self, path):
        self.path = path
        self.is_geodesic = self.is_path_too_short()

        if self.is_geodesic:
            self.joints = []
            return

        self.joints = [
            Joint.from_vertices(self.tri,
                                self.path[i-1],
                                self.path[i],
                                self.path[i+1],
                                i)
            for i in range(1, len(path) - 1)
        ]
        self.joints.extend([j.twin for j in self.joints])

    def update_path(self, old_joint, bypass):
        # remove b
        b_index = old_joint.path_data
        del self.path[b_index]
        if b_index < 0:
            b_index += len(self.path) + 1

        # insert bypass to replace it
        if not old_joint.direction == Joint.Direction.FORTH:
            bypass.reverse()

        self.path[b_index:b_index] = bypass[1:-1]

        # could be made much more efficient, but for that we'd need to compute the indices correctly, and
        # have every joint keep pointers to its next & previous joints
        self.set_path(self.path)

    def get_minimal_wedge_angle(self):
        if self.is_geodesic:
            return np.inf
        return self.get_minimal_joint().wedge_angle

    def get_minimal_joint(self):
        return min(self.joints, key=lambda j: j.wedge_angle)

    def flipout_the_minimal_wedge(self):
        if self.is_geodesic:
            return

        joint = self.get_minimal_joint()
        if is_reflex_or_flat(joint.wedge_angle, self.reflex_angle_threshold_for_flip_out):
            self.is_geodesic = True
            return

        bypass = self.flipout(joint)
        self.update_path(joint, bypass)

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

        if is_reflex_or_flat(min_angle, threshold=self.reflex_angle_threshold_for_flip_out):
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
        if len(self.path) == 0:
            return []
        return self.path + [self.path[0]]

    @property
    def length(self):
        return np.sum([self.tri.get_edge(self.path[i - 1], self.path[i]).length for i in range(len(self.path))])

    def is_path_too_short(self):
        return len(self.path) <= 1

    def set_path(self, path):
        if len(path) > 0 and path[0] == path[-1]:
            path.pop()
        self.path = path
        self.is_geodesic = self.is_path_too_short()

        if self.is_geodesic:
            self.joints = []
            return

        self.joints = [
            Joint.from_vertices(self.tri,
                                self.path[i - 2],
                                self.path[i - 1],
                                self.path[i],
                                i - 1)
            for i in range(len(path))
        ]
        self.joints.extend([j.twin for j in self.joints])

    def flipout_the_minimal_wedge(self):
        super().flipout_the_minimal_wedge()
        if len(self.path) == 2 and len(self.tri.get_all_edges_between(self.path[0], self.path[1])) == 1:
            self.path.pop()


class MultiplePathShortener(BaseShortener):
    def __init__(self, triangulation):
        super().__init__(triangulation)
        self.shorteners = []
        self.fixed_vertices = []

    def get_path(self):
        return [shortener.get_path() for shortener in self.shorteners]

    @staticmethod
    def split_paths(paths):
        while MultiplePathShortener.split_paths_one_at_a_time(paths):
            pass

    @staticmethod
    def split_paths_one_at_a_time(paths):
        for i in range(len(paths)):
            path = paths[i]
            for other_path in paths:
                if path == other_path:
                    continue

                for pos in range(1, len(path) - 1):
                    if path[pos] in other_path:
                        path1 = path[:pos+1]
                        path2 = path[pos:]
                        del paths[i]
                        if path1[0] == path2[-1]:
                            # path was a loop
                            path2.pop()
                            path2.extend(path1)
                            path2.append(None)
                            paths.append(path2)
                        else: # path
                            paths.append(path1)
                            paths.append(path2)
                        return True

        return False

    def set_path(self, paths):
        # split paths that contain points from other paths
        self.split_paths(paths)

        # create lots of shorteners
        for path in paths:
            if path[-1] is None:
                roi = ROI.PATH
                path.pop()
            else:
                roi = ROI.LOOP if path[0] == path[-1] else ROI.PATH
            shortener = get_shortener(roi, self.tri)
            shortener.set_path(path)
            self.shorteners.append(shortener)

    def flipout_the_minimal_wedge(self):
        if self.is_geodesic:
            return

        shortener = min(self.shorteners, key=lambda shortener: shortener.get_minimal_wedge_angle())
        if shortener.is_geodesic:
            self.is_geodesic = True
            return

        shortener.flipout_the_minimal_wedge()

