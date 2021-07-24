import pyvista as pv
import numpy as np
from tkinter import *
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
        self.mesh_obj = self.set_up_extrinsic_mesh(url)
        self.tri = self.set_up_triangulation()
        self.remove_extrinsic_mesh_faces()
        self.tri_obj = self.set_up_intrinsic_triangulation()
        self.plotter = self.set_up_plotter()
        self.path_shortener = PathShortener(self.tri)
        self.path_picker = self.set_up_path_picker()
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

    def set_up_extrinsic_mesh(self, url):
        obj = pv.read(url)

        if not obj.is_all_triangles():
            obj.triangulate()

        self.mesh_obj = obj
        return obj

    def set_up_triangulation(self):
        V = self.mesh_obj.points
        F = self.mesh_obj.faces
        F = np.reshape(F, (len(F) // 4, 4))
        F = np.delete(F, 0, axis=1)
        tri = Triangulation(V, F)
        tri.init_coloring(len(F), prefer.TRIANGULATION_FACES_COLOR_MAP_SIZE)

        self.tri = tri
        return tri

    def set_up_intrinsic_triangulation(self):
        V, F = self.tri.get_faces()
        obj = pv.PolyData(V, F)

        obj.cell_arrays['TriagColoring'] = self.tri.get_coloring()

        self.tri_obj = obj
        return obj

    def remove_extrinsic_mesh_faces(self):
        V, E = self.mesh_obj.points, self.mesh_obj.faces
        self.mesh_obj = pv.PolyData(V, lines=E)

    def set_up_plotter(self):
        self.plotter = pv.Plotter()

        self.add_extrinsic_mesh()
        self.add_intrinsic_triangulation()

        return self.plotter

    def add_extrinsic_mesh(self):
        mesh_kwargs = {
            'name': 'mesh_edges',
            'color': prefer.MESH_EDGE_COLOR,
            'line_width': prefer.MESH_EDGE_WIDTH,
            'show_edges': prefer.SHOW_MESH_EDGES,
        }
        self.plotter.add_mesh(self.mesh_obj, **mesh_kwargs)

    def add_intrinsic_triangulation(self):
        self.set_up_intrinsic_triangulation()
        tri_kwargs = {
            'edge_color': prefer.TRIANGULATION_EDGE_COLOR,
            'line_width': prefer.TRIANGULATION_EDGE_WIDTH,
            'show_edges': prefer.SHOW_TRIANGULATION_EDGES,
            'show_scalar_bar': False,
            'render_points_as_spheres': True,
        }
        if prefer.SHOW_TRIANGULATION_FACES:
            tri_kwargs['cmap'] = prefer.TRIANGULATION_FACES_COLOR_MAP
            self.tri_obj.set_active_scalars('TriagColoring')
        else:
            tri_kwargs['color'] = 'Gold'
            self.tri_obj.set_active_scalars(None)

        self.plotter.remove_actor('tri')
        self.plotter.add_mesh(self.tri_obj, **tri_kwargs, name='tri')

    def show(self):
        self.plotter.show()

    def set_up_path_picker(self):
        path_picker = PathPicker(self, line_width=prefer.PICKED_PATH_WIDTH, point_size=prefer.PICKED_POINT_SIZE)

        self.path_picker = path_picker
        return path_picker

    def set_up_events(self):
        bindings = {
            prefer.KEY_EVENT_SINGLE_SOURCE_DIJKSTRA: self.on_single_source,
            prefer.KEY_EVENT_MAKE_GEODESIC: self.on_make_geodesic,
            prefer.KEY_EVENT_FLIPOUT: self.on_flip_out,
            prefer.KEY_EVENT_EDGE_FLIP: self.on_flip_edge,
            prefer.KEY_EVENT_SHOW_INFO: self.on_info,
            prefer.KEY_EVENT_CLEAR_PICKED_PATH: self.path_picker.on_clear,
            prefer.KEY_EVENT_UNDO_PICK: self.path_picker.on_undo,
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
        self.path_shortener.flipout_the_minimal_wedge()

    def on_flip_edge(self):
        e = self.path_picker.get_corresponding_edge(self.tri)
        if e is None:
            return

        self.path_picker.on_clear()

        self.tri.flip(e)
        self.add_intrinsic_triangulation()

    def on_info(self):
        pass


scene = Scene('C:/Users/Arnon/Desktop/cup3.obj')
scene.show()





exit()

from pyvista import examples

# Load a global topography surface and decimate it
land = examples.download_topo_land().triangulate().decimate(0.98)

cape_town = land.find_closest_point((0.790801, 0.264598, -0.551942))
dubai = land.find_closest_point((0.512642, 0.745898, 0.425255))
bangkok = land.find_closest_point((-0.177077, 0.955419, 0.236273))
rome = land.find_closest_point((0.718047, 0.163038, 0.676684))

a = land.geodesic(cape_town, dubai)
b = land.geodesic(cape_town, bangkok)
c = land.geodesic(cape_town, rome)



def picker(event):
    print("A")

# mesh points
vertices = np.array([[0, 0, 0],
                     [1, 0, 0],
                     [1, 1, 0],
                     [0, 1, 0],
                     [0.5, 0.5, -1]])

# mesh faces
faces = np.hstack([[4, 0, 1, 2, 3],  # square
                   [3, 0, 1, 4],     # triangle
                   [3, 1, 2, 4]])    # triangle

surf = pv.PolyData(vertices, faces)

plotter = pv.Plotter()
# plotter.add_mesh(surf)
plotter.add_mesh(a+b+c, line_width=10, color="red", label="Geodesic Path")
plotter.add_mesh(land, show_edges=True)
plotter.add_legend()

plotter.enable_path_picking(callback=picker, show_message=False)
plotter.show()

# plot each face with a different color
# surf.plot(scalars=np.arange(3), cpos=[-1, 1, 0.5])

