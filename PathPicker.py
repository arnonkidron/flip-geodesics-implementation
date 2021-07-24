import pyvista as pv
import ViewPreferences as prefer


class PathPicker:
    def __init__(self, scene, tolerance=0.025, **kwargs):
        self.path_actor = pv.PolyData()
        self.indices = []
        self.whole_path_indices = []

        self.scene = scene
        self.show_path = True
        self.show_endpoints = True
        self.show_midpoints = True
        self.allow_pick_intersection_points_of_mesh_and_triangulation = False

        self.kwargs = kwargs
        self.kwargs.setdefault('color', prefer.PICKED_PATH_COLOR)
        self.kwargs.setdefault('point_size', prefer.PICKED_POINT_SIZE)
        self.kwargs.setdefault('line_width', prefer.PICKED_PATH_WIDTH)

        self.kwargs.setdefault('render_points_as_spheres', True)
        self.kwargs.setdefault('pickable', False)
        self.kwargs.setdefault('reset_camera', False)
        self.path_name = '_picked_path'

        self.scene.plotter.enable_point_picking(
            callback=self.on_pick,
            use_mesh=False,
            tolerance=tolerance,
            show_point=False,
            show_message=False,
        )

    def is_empty(self):
        return len(self.indices) == 0

    @property
    def last_index(self):
        return self.indices[-1]

    def get_single_point_index(self):
        if len(self.indices) != 1:
            return None

        return self.indices[0]

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
        self.pick_edge(e.next)

    def on_pick_twin_edge(self):
        e = self.get_corresponding_edge()
        if e is None:
            return
        self.pick_edge(e.twin)

    def pick_edge(self, e):
        self.indices = [e.origin, e.dst]
        self.reconstruct_by_indices()
        self.scene.on_info()

    def is_intersection_point(self, idx):
        return idx >= self.scene.mesh_actor.n_points

    def on_pick(self, input_point):
        point_finder = self.scene.mesh_actor
        if self.allow_pick_intersection_points_of_mesh_and_triangulation:
            point_finder = self.scene.tri_actor
        idx = point_finder.find_closest_point(input_point)
        point = point_finder.points[idx]

        if self.allow_pick_intersection_points_of_mesh_and_triangulation:
            # if the current or the last point is such,
            # then clear the path before moving on
            if self.is_intersection_point(idx) \
                    or (not self.is_empty() and self.is_intersection_point(self.last_index)):
                self.on_clear()

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

    def init_path_vertices(self):
        last_point = self.scene.tri_actor.points[self.last_index]
        first_point = self.scene.tri_actor.points[self.indices[0]]
        self.path_actor = pv.PolyData([last_point, first_point])
        self.whole_path_indices = [self.indices[0]]

    def update_path_edges(self, prev_idx, current_idx):
        new_part = self.scene.tri.find_shortest_path(prev_idx, current_idx)
        self.whole_path_indices.pop()
        self.whole_path_indices.extend(new_part)

        V = self.scene.tri.V[new_part]
        E = [2 * len(V) - 1]
        rang = list(range(len(V)))
        E += rang
        rang.reverse()
        E += rang
        poly_data = pv.PolyData(V, lines=E)
        self.path_actor += poly_data

    def update_path_vertices(self, current_point):
        if self.show_midpoints:
            self.path_actor += pv.PolyData(current_point)
        elif self.show_endpoints:
            self.path_actor.points[0] = current_point
        elif len(self.indices) == 2:
            self.path_actor = pv.PolyData()

    def add_actor(self):
        if self.show_path:
            self.remove_actor()
            self.scene.plotter.add_mesh(self.path_actor, **self.kwargs, name=self.path_name)

    def remove_actor(self):
        self.scene.plotter.remove_actor(self.path_name)

    def on_clear(self):
        self.path_actor = pv.PolyData()
        self.indices = []
        self.whole_path_indices = []
        self.remove_actor()
        self.scene.remove_text()
        return

    def on_undo(self):
        if self.is_empty():
            return

        self.indices.pop()
        if self.is_empty():
            self.on_clear()
            return

        self.reconstruct_by_indices()
        self.scene.remove_text()

    def reconstruct_by_indices(self):
        self.path_actor = pv.PolyData()
        if self.show_endpoints or len(self.indices) == 1:
            self.init_path_vertices()

        for i in range(len(self.indices) - 1):
            self.update_path_edges(self.indices[i], self.indices[i + 1])
            if self.show_midpoints:
                current_point = self.scene.tri_actor.points[self.indices[i]]
                self.path_actor += pv.PolyData(current_point)

        self.add_actor()








