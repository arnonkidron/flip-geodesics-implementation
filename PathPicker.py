import pyvista as pv
import ViewPreferences as prefer
import exceptions
import numpy as np
from utils import ROI


class PathVisualizer:
    def __init__(self, scene, **kwargs):
        self.path_actor = None
        self.fixed_points_actor = None
        self.indices = []
        self.whole_path_indices = []
        self.fixed_positions = []

        self.roi = None

        self.scene = scene

        self.kwargs = kwargs
        self.kwargs.setdefault('color', prefer.RESULT_PATH_COLOR)
        self.kwargs.setdefault('point_size', prefer.PATH_POINT_SIZE)
        self.kwargs.setdefault('line_width', prefer.PATH_EDGE_WIDTH)

        self.kwargs.setdefault('render_points_as_spheres', True)
        self.kwargs.setdefault('render_lines_as_tubes', True)
        self.kwargs.setdefault('pickable', False)
        self.kwargs.setdefault('reset_camera', False)
        self.kwargs.setdefault('name', '_result_path')

        self.show_path = prefer.SHOW_RESULT_PATH
        self.show_path_end_points = prefer.SHOW_RESULT_PATH_END_POINTS
        self.show_path_all_points = prefer.SHOW_RESULT_PATH_ALL_POINTS
        self.show_network_fixed_points = prefer.SHOW_NETWORK_FIXED_POINTS

    def is_empty(self):
        if self.indices is None:
            return False
        return len(self.indices) == 0

    @property
    def first_index(self):
        return self.indices[0]

    @property
    def last_index(self):
        return self.indices[-1]

    def get_path(self):
        return self.indices

    def get_path_edge_tuples_set(self):
        if self.is_empty():
            return set()

        tuples = [(self.whole_path_indices[i], self.whole_path_indices[i+1]) for i in range(len(self.whole_path_indices) - 1)]
        if self.roi == ROI.LOOP:
            tuples.append((self.whole_path_indices[-1], self.whole_path_indices[0]))

        tuples.extend([(v, u) for (u, v) in tuples])
        return set(tuples)

    def set_path(self, path):
        self.indices = path

    def set_path_as_one_edge(self, e):
        self.set_path([e.origin, e.dst])
        self.roi = ROI.PATH

    def add_actor(self):
        if not self.show_path or self.is_empty():
            self.remove_actor()
            return

        self.reconstruct_by_indices()

        self.scene.plotter.add_mesh(self.path_actor, **self.kwargs)

    def remove_actor(self):
        self.scene.plotter.remove_actor(self.kwargs['name'])

    def clear(self):
        self.path_actor = pv.PolyData()
        self.indices = []
        self.whole_path_indices = []
        self.remove_actor()

    def init_path_vertices(self):
        last_point = self.scene.V[self.last_index]
        first_point = self.scene.V[self.first_index]
        self.path_actor = pv.PolyData([last_point, first_point])
        self.fixed_points_actor = pv.PolyData()
        self.whole_path_indices = [self.first_index]
        self.roi = ROI.VERTEX

    def update_path_edges(self, prev_idx, current_idx):
        new_part = self.scene.tri.find_shortest_path(prev_idx, current_idx)
        self.whole_path_indices.pop()
        common_vertices = list(set.intersection(set(self.whole_path_indices), set(new_part)))
        prev_len = len(self.whole_path_indices)
        self.whole_path_indices.extend(new_part)

        if len(common_vertices) == 0:
            self.roi = ROI.PATH
        elif len(common_vertices) == 1 and common_vertices[0] == self.first_index and self.roi == ROI.PATH:
            self.roi = ROI.LOOP
        else:
            self.roi = ROI.NETWORK
            for new_index in common_vertices:
                if new_index in self.whole_path_indices[:prev_len]:
                    self.fixed_positions.append(new_index)
                    self.fixed_points_actor += pv.PolyData(self.scene.tri.V[new_index])

        # add to actor
        V = self.scene.tri.V[new_part]
        E = [None]
        last_vertex = len(V) - 1
        for i in range(last_vertex):
            E.append(i)
            e = self.scene.tri.get_edge(new_part[i], new_part[i+1])
            if e is None:
                raise exceptions.NonExistentEdgeException(new_part[i], new_part[i+1])

            if self.scene.slow_edge is not None:
                continue

            if self.roi == ROI.LOOP and len(self.indices) == 3:
                # parallel edges, now we must choose the other one
                e = self.scene.tri.get_all_edges_between(new_part[i], new_part[i+1])[-1]

            intersections = e.get_intersections(self.scene.tri.mesh)
            if e.intersections_status != e.Status.FAILED \
                    and len(intersections) > 0:
                index_begin = len(V)
                V = np.vstack((V, [p.coords for p in intersections]))
                index_end = len(V)
                E.extend(list(range(index_begin, index_end)))
            elif e.intersections_status == e.Status.FAILED:
                if not prefer.SHOW_FAILED_PATH_EDGES:
                    E[0] = len(E) - 1
                    poly_data = pv.PolyData(V, lines=E)
                    self.path_actor += poly_data
                    E = [None]
                    continue
                index_begin = len(V)
                intersection_coords = [p.coords for p in intersections]
                twin_intersection_coords = [p.coords for p in e.twin.get_intersections(self.scene.tri.mesh)]
                V = np.vstack((V, intersection_coords))
                if len(twin_intersection_coords) > 0:
                    V = np.vstack((V, twin_intersection_coords[::-1]))
                index_end = len(V)
                E.extend(list(range(index_begin, index_end)))

        E.append(last_vertex)

        E[0] = len(E) - 1
        poly_data = pv.PolyData(V, lines=E)
        self.path_actor += poly_data

    def update_path_vertices(self, current_point):
        if self.show_path_all_points:
            self.path_actor += pv.PolyData(current_point)
        elif self.show_path_end_points:
            self.path_actor.points[0] = current_point
        elif len(self.indices) == 2:
            self.path_actor = pv.PolyData()

    def reconstruct_by_indices(self):
        if self.is_empty():
            return

        self.path_actor = pv.PolyData()
        if self.show_path_end_points or len(self.indices) == 1:
            self.init_path_vertices()

        for i in range(len(self.indices) - 1):
            self.update_path_edges(self.indices[i], self.indices[i + 1])
            if self.show_path_all_points:
                current_point = self.scene.V[self.indices[i]]
                self.path_actor += pv.PolyData(current_point)


class PathPicker(PathVisualizer):
    def __init__(self, scene, tolerance=0.025, **kwargs):
        kwargs.setdefault('color', prefer.PICKED_PATH_COLOR)
        kwargs.setdefault('name', '_picked_path')

        super().__init__(scene, **kwargs)

        self.show_path = prefer.SHOW_PICKED_PATH
        self.show_path_end_points = prefer.SHOW_PICKED_PATH_END_POINTS
        self.show_path_all_points = prefer.SHOW_PICKED_PATH_ALL_POINTS

        self.scene.plotter.enable_point_picking(
            callback=self.on_pick,
            use_mesh=False,
            tolerance=tolerance,
            show_point=False,
            show_message=False,
        )

    def add_actor(self):
        if not self.show_path:
            self.remove_actor()
            return

        kwargs = self.kwargs.copy()
        if len(self.indices) == 1 and self.is_intersection_point(self.last_index):
                kwargs['color'] = prefer.PICKED_INTERSECTION_POINTS_COLOR
        self.scene.plotter.add_mesh(self.path_actor, **kwargs)

        if self.fixed_points_actor.n_points > 0:
            self.scene.plotter.add_mesh(self.fixed_points_actor,
                                    name='_fixed_points', render_points_as_spheres=True,
                                    color=prefer.NETWORK_FIXED_POINTS_COLOR, point_size=self.kwargs['point_size'])

    def set_path(self, path):
        self.indices = path
        self.reconstruct_by_indices()
        self.add_actor()

    def get_single_point_index(self):
        if len(self.indices) != 1:
            return None

        return self.first_index

    def get_corresponding_edge(self):
        """
        if the picked path is an edge in the triangulation, return it
        else, None
        """
        if len(self.indices) != 2:
            return None

        return self.scene.tri.get_edge(self.indices[0], self.indices[1])

    def on_pick_next_edge(self):
        e = self.get_corresponding_edge()
        if e is None:
            return
        self.set_path_as_one_edge(e.next)
        self.scene.on_info()

    def on_pick_twin_edge(self):
        e = self.get_corresponding_edge()
        if e is None:
            return
        self.set_path_as_one_edge(e.twin)
        self.scene.on_info()

    def is_intersection_point(self, idx):
        return idx >= self.scene.mesh_actor.n_points

    def on_pick(self, input_point):
        point_finder = self.scene.mesh_actor
        if prefer.ALLOW_PICK_INTERSECTION_POINTS:
            point_finder = self.scene.tri_edge_actor
        idx = point_finder.find_closest_point(input_point)
        point = point_finder.points[idx]

        if prefer.ALLOW_PICK_INTERSECTION_POINTS:
            # if the current or the last point is such,
            # then clear the path before moving on
            if self.is_intersection_point(idx) \
                    or (not self.is_empty() and self.is_intersection_point(self.last_index)):
                self.clear()

        if self.is_empty():
            self.indices = [idx]
            self.init_path_vertices()
        else:
            if idx == self.last_index:
                return
            self.update_path_vertices(point)
            self.update_path_edges(self.last_index, idx)
            self.indices.append(idx)

        self.add_actor()

    def on_undo(self):
        if self.is_empty():
            return

        self.indices.pop()
        if self.is_empty():
            self.clear()
            return

        self.reconstruct_by_indices()
        self.add_actor()
        self.scene.remove_text()

    def on_close_loop(self):
        self.on_pick(self.scene.tri.V[self.first_index])
        self.is_loop = True

