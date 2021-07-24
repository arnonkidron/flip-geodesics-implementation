import pyvista as pv
import numpy as np


class PathPicker:
    def __init__(self, scene, tolerance=0.025, **kwargs):
        self.path_obj = pv.PolyData()
        self.indices = []

        self.scene = scene
        self.show_path = True
        self.show_endpoints = True
        self.show_midpoints = False

        self.kwargs = kwargs
        self.kwargs.setdefault('pickable', False)
        self.kwargs.setdefault('color', 'pink')
        self.kwargs.setdefault('point_size', 50)
        # self.kwargs.setdefault('line_width', 5)

        self.kwargs.setdefault('reset_camera', False)
        self.path_name = '_picked_path'

        self.scene.plotter.add_key_event('c', self.on_clear)
        self.scene.plotter.add_key_event('z', self.on_undo)
        self.scene.plotter.enable_point_picking(callback=self.on_pick, use_mesh=False,
                                          tolerance=tolerance, show_point=False)

    def is_empty(self):
        return len(self.indices) == 0

    @property
    def last_index(self):
        return self.indices[-1]

    def get_single_point_index(self):
        if len(self.indices) != 1:
            return None

        return self.indices[0]

    def get_corresponding_edge(self, triangulation):
        """
        Input: a triangulation object
        Output: if then picked path is an edge in the triangulation, return it
                else, None
        """
        if len(self.indices) != 2:
            return None

        return triangulation.get_edge(self.indices[0], self.indices[1])

    def on_pick(self, input_point):
        idx = self.scene.mesh_obj.find_closest_point(input_point)

        point = self.scene.tri_obj.points[idx]
        if self.is_empty():
            self.indices = [idx]
            self.init_path_vertices()
        else:
            self.update_path_vertices(point)
            self.update_path_edges(self.last_index, idx)
            self.indices.append(idx)

        self.show()

    def init_path_vertices(self):
        last_point = self.scene.tri_obj.points[self.last_index]
        first_point = self.scene.tri_obj.points[self.indices[0]]
        self.path_obj = pv.PolyData([last_point, first_point])

    def update_path_edges(self, prev_idx, current_idx):
        V = self.scene.tri.find_shortest_path(prev_idx, current_idx)
        E = [len(V)] + list(range(len(V)))
        poly_data = pv.PolyData(V, lines=E)
        self.path_obj += poly_data

    def update_path_vertices(self, current_point):
        if self.show_midpoints:
            self.path_obj += pv.PolyData(current_point)
        elif self.show_endpoints:
            self.path_obj.points[0] = current_point
        elif len(self.indices) == 2:
            self.path_obj = pv.PolyData()

    def show(self):
        if self.show_path:
            self.hide()
            self.scene.plotter.add_mesh(self.path_obj, **self.kwargs, name=self.path_name)

    def hide(self):
        self.scene.plotter.remove_actor(self.path_name)

    def on_clear(self):
        self.path_obj = pv.PolyData()
        self.indices = []
        self.hide()
        return

    def on_undo(self):
        if self.is_empty():
            return

        self.indices.pop()
        if self.is_empty():
            self.hide()
            return

        self.path_obj = pv.PolyData()
        if self.show_endpoints or len(self.indices) == 1:
            self.init_path_vertices()

        for i in range(len(self.indices) - 1):
            self.update_path_edges(self.indices[i], self.indices[i + 1])

        self.show()








