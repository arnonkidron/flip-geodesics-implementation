import pyvista as pv
import numpy as np
from tkinter import Tk, messagebox
from tkinter.filedialog import askopenfilename
from FlipEdgeNetwork import *
from Triangulation import *
from edge import *
from Solver import *
from PathPicker import *
from PathShortener import *
import ViewPreferences as prefer


class Scene:
    def __init__(self, url=None):
        if url is None:
            url = self.ask_for_url()

        self.plotter = pv.Plotter()
        self.mesh_actor = self.set_up_extrinsic_mesh(url)
        self.tri = self.set_up_triangulation()
        self.tri_actor = None

        self.slow_generator = None

        self.add_extrinsic_mesh_actor()
        self.add_intrinsic_triangulation()

        self.path_shortener = PathShortener(self.tri)
        self.path_picker = self.set_up_path_picker()
        self.text_actor = None
        self.set_up_events()

    @staticmethod
    def ask_for_url():
        # we don't want a full GUI, so keep the root window from appearing
        Tk().withdraw()

        # show an "Open" dialog box and return the path to the selected file
        mesh_url = askopenfilename(
            title='Choose a mesh file to load',
            filetypes=[
                ('All Files (*.*)', '*'),
            ]
        )

        if mesh_url == '':
            exit("flip_geodesic load error: No file chosen")

        return mesh_url

    @staticmethod
    def warn(title, msg):
        """
        Outcome: display the message in a tkinter window, with retry & cancel
        Output: whether the user has chosen to retry
        """
        Tk().withdraw()
        return messagebox.showwarning(title=title, message=msg, type=messagebox.OK)

    def show(self):
        self.plotter.show()

    def set_up_extrinsic_mesh(self, url):
        actor = pv.read(url)

        if hasattr(actor, 'is_all_triangles'):
            if not actor.is_all_triangles():
                actor.triangulate()

        V = actor.points
        if hasattr(actor, 'faces'):
            E = actor.faces
        else:
            E = actor.cells
        self.mesh_actor = pv.PolyData(V, lines=E)
        return self.mesh_actor

    def set_up_triangulation(self):
        V = self.mesh_actor.points
        F = self.mesh_actor.lines
        F = np.reshape(F, (len(F) // 4, 4))
        F = np.delete(F, 0, axis=1)
        tri = Triangulation(V, F)
        tri.init_coloring(len(F), prefer.TRIANGULATION_FACES_COLOR_MAP_SIZE)

        self.tri = tri
        return tri

    def add_extrinsic_mesh_actor(self):
        mesh_kwargs = {
            'name': 'mesh_edges',
            'color': prefer.MESH_EDGE_COLOR,
            'line_width': prefer.MESH_EDGE_WIDTH,
            'show_edges': prefer.SHOW_MESH_EDGES,
        }
        self.plotter.add_mesh(self.mesh_actor, **mesh_kwargs)

    def add_intrinsic_triangulation(self):
        # advance intersection point computation by one step
        if self.slow_generator is not None:
            self.plotter.remove_actor('slow_edge')
            self.plotter.remove_actor('hit_edge')
            self.plotter.remove_actor('hit_point')
            self.plotter.remove_actor('next_vec')

            point, mesh_edge, tri_edge, new_vec = next(self.slow_generator)
            if point is None:
                self.slow_generator = None
            else:
                actor = pv.PolyData([self.mesh_actor.points[tri_edge.origin], self.mesh_actor.points[tri_edge.dst]], [2, 0,1])
                # self.plotter.add_mesh(actor, name='slow_edge', style='wireframe', color=prefer.TRIANGULATION_EDGE_COLOR, render_lines_as_tubes=True, line_width=prefer.PICKED_PATH_WIDTH)
                actor = pv.PolyData([self.mesh_actor.points[mesh_edge.origin], self.mesh_actor.points[mesh_edge.dst]], [2, 0,1])
                self.plotter.add_mesh(actor, name='hit_edge', style='wireframe', color=prefer.COMPUTED_INTERSECTING_EDGE_COLOR, render_lines_as_tubes=True, line_width=prefer.PICKED_PATH_WIDTH)
                self.plotter.add_mesh(pv.PolyData(point), name='hit_point', color=prefer.COMPUTED_INTERSECTION_POINTS_COLOR, render_points_as_spheres=True, point_size=prefer.PICKED_POINT_SIZE)

                actor = pv.PolyData([point, point + new_vec], [2, 0,1])
                self.plotter.add_mesh(actor, name='next_vec', style='wireframe',
                                      color='Green', render_lines_as_tubes=True,
                                      line_width=prefer.PICKED_PATH_WIDTH)

        # set up
        V, F, coloring = self.tri.get_poly_data()
        self.tri_actor = pv.PolyData(V, F)
        self.tri_actor.cell_arrays['TriagColoring'] = coloring

        tri_kwargs = {
            'edge_color': prefer.TRIANGULATION_EDGE_COLOR,
            'line_width': prefer.TRIANGULATION_EDGE_WIDTH,
            'show_edges': prefer.SHOW_TRIANGULATION_EDGES,
            'show_scalar_bar': False,
            'render_points_as_spheres': True,
        }
        if prefer.SHOW_TRIANGULATION_FACES:
            tri_kwargs['cmap'] = prefer.TRIANGULATION_FACES_COLOR_MAP
            self.tri_actor.set_active_scalars('TriagColoring')
        else:
            tri_kwargs['color'] = 'Gold'
            self.tri_actor.set_active_scalars(None)

        self.plotter.remove_actor('tri')
        self.plotter.add_mesh(self.tri_actor, **tri_kwargs, name='tri')

    def set_up_path_picker(self):
        path_picker = PathPicker(self, line_width=prefer.PICKED_PATH_WIDTH, point_size=prefer.PICKED_POINT_SIZE)

        self.path_picker = path_picker
        return path_picker

    def set_up_events(self):
        bindings = {
            prefer.KEY_EVENT_SINGLE_SOURCE_DIJKSTRA: self.on_single_source,
            prefer.KEY_EVENT_MAKE_GEODESIC: self.on_make_geodesic,
            prefer.KEY_EVENT_FLIPOUT: self.on_flip_out,
            prefer.KEY_EVENT_EDGE_FLIP: self.on_edge_flip,
            prefer.KEY_EVENT_SHOW_INFO: self.on_info,
            prefer.KEY_EVENT_CLEAR_PICKED_PATH: self.path_picker.on_clear,
            prefer.KEY_EVENT_UNDO_PICK: self.path_picker.on_undo,
            prefer.KEY_EVENT_PICK_NEXT_EDGE: self.path_picker.on_pick_next_edge,
            prefer.KEY_EVENT_PICK_TWIN_EDGE: self.path_picker.on_pick_twin_edge,
            prefer.KEY_EVENT_RE_RENDER: self.add_intrinsic_triangulation,
        }
        for key, callback in bindings.items():
            self.plotter.add_key_event(key, callback)

    def on_single_source(self):
        idx = self.path_picker.get_single_point_index()
        if idx is None:
            return

        pass

    def on_make_geodesic(self):
        path = self.path_picker.whole_path_indices
        self.path_shortener.set_path(path)
        self.path_shortener.make_geodesic()

    def on_flip_out(self):
        path = self.path_picker.whole_path_indices
        self.path_shortener.set_path(path)
        tmp = self.path_shortener.flipout_the_minimal_wedge()

        tmp_picker = PathPicker(self, color='Red')
        tmp_picker.path_name = "_bypass_path"
        for v in tmp:
            tmp_picker.on_pick(self.tri.V[v])



    def on_edge_flip(self):
        old_edge = self.path_picker.get_corresponding_edge()
        if old_edge is None:
            return

        self.path_picker.on_clear()

        new_edge = self.tri.flip(old_edge)

        # check for error
        if isinstance(new_edge, str):
            self.warn(title="Edge flip fail", msg=new_edge)
            return

        # compute intersection points
        if prefer.COMPUTE_INTERSECTION_POINTS_ONE_AT_A_TIME:
            self.slow_generator = new_edge.init_intersections_one_at_a_time()
        else:
            new_edge.init_intersections()

        self.add_intrinsic_triangulation()

    def on_info(self):
        msg = None

        idx = self.path_picker.get_single_point_index()
        if idx is not None:
            coords = self.tri_actor.points[idx]
            if self.path_picker.is_intersection_point(idx):
                msg = "Intersection point\n" \
                      "({:.4f},{:.4f},{:.4f})\n"\
                    .format(coords[0], coords[1], coords[2])
            else:
                deg = len(self.tri.in_edges[idx])
                msg = "Vertex #{}\n" \
                      "({:.4f},{:.4f},{:.4f})\n" \
                      "Degree {}\n"\
                    .format(idx, coords[0], coords[1], coords[2], deg)

        e = self.path_picker.get_corresponding_edge()
        if e is not None:
            msg = e.get_info()

        self.remove_text()
        if msg is None:
            self.text_actor = None
        else:
            self.text_actor = self.plotter.add_text(msg, name='info')
            print(msg)

    def remove_text(self):
        self.plotter.remove_actor(self.text_actor)


if __name__ == '__main__':
    scene = Scene()
    # scene = Scene()
    scene.show()

