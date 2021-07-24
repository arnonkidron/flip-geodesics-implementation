import pyvista as pv
import numpy as np
from tkinter import *
from tkinter.filedialog import askopenfilename
from FlipEdgeNetwork import *
from Triangulation import *
from edge import *
from Solver import *
from PathPicker import *

############################
#    view preferences
############################
MESH_EDGE_WIDTH = 0.5
TRIANGULATION_EDGE_WIDTH = 4
PICKED_PATH_WIDTH = 30
SHOW_MESH_EDGES = True
SHOW_TRIANGULATION_EDGES = False
SHOW_TRIANGULATION_FACES = False

TRIANGULATION_FACES_COLOR_MAP = 'Accent'
TRIANGULATION_FACES_COLOR_MAP_SIZE = 8


class Scene:
    def __init__(self, url=None):
        if url is None:
            url = self.ask_for_url()
        self.mesh_obj = self.set_up_extrinsic_mesh(url)
        self.tri = self.set_up_triangulation()
        self.tri_obj = self.set_up_intrinsic_mesh()
        self.plotter = self.set_up_plotter()
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
        mesh = pv.read(url)

        if not mesh.is_all_triangles():
            mesh.triangulate()

        self.mesh_obj = mesh
        return mesh

    def set_up_triangulation(self):
        V = self.mesh_obj.points
        F = self.mesh_obj.faces
        F = np.reshape(F, (len(F) // 4, 4))
        F = np.delete(F, 0, axis=1)
        tri = Triangulation(V, F)
        self.tri = tri
        return tri

    def set_up_intrinsic_mesh(self):
        V, F = self.tri.get_faces()
        obj = pv.PolyData(V, F)

        obj.cell_arrays['TriagColor'] = self.tri.get_coloring(TRIANGULATION_FACES_COLOR_MAP_SIZE)

        self.tri_obj = obj
        return obj

    def set_up_plotter(self):
        plotter = pv.Plotter()

        tri_kwargs = {
            'name': 'tri',
            'edge_color': 'Black',
            'line_width': TRIANGULATION_EDGE_WIDTH,
            'show_edges': SHOW_TRIANGULATION_EDGES,
            'show_scalar_bar': False,
        }
        mesh_kwargs = {
            'name': 'mesh',
            'edge_color': 'Black',
            'line_width': MESH_EDGE_WIDTH,
            'show_edges': SHOW_MESH_EDGES,
            'show_scalar_bar': False,
        }
        if SHOW_TRIANGULATION_FACES:
            if SHOW_MESH_EDGES:
                mesh_kwargs['style'] = 'wireframe'
                mesh_kwargs['color'] = 'Black'
            tri_kwargs['cmap'] = TRIANGULATION_FACES_COLOR_MAP
            self.tri_obj.set_active_scalars('TriagColor')
        else:
            mesh_kwargs['color'] = 'Gold'
            tri_kwargs['color'] = 'Gold'
            if SHOW_TRIANGULATION_EDGES:
                tri_kwargs['color'] = 'Black'
                tri_kwargs['style'] = 'wireframe'
            self.tri_obj.set_active_scalars(None)

        # plotter.add_mesh(self.mesh_obj, **mesh_kwargs)
        plotter.add_mesh(self.tri_obj, **tri_kwargs)

        self.plotter = plotter
        return plotter

    def show(self):
        self.plotter.show()

    def set_up_path_picker(self):
        path_picker = PathPicker(self, line_width=PICKED_PATH_WIDTH)

        self.path_picker = path_picker
        return path_picker

    def set_up_events(self):
        self.plotter.add_key_event('[', self.on_single_source)
        self.plotter.add_key_event('o', self.on_make_geodesic)
        self.plotter.add_key_event('i', self.on_flip_out)
        self.plotter.add_key_event('u', self.on_flip_edge)

    def on_single_source(self):
        idx = self.path_picker.get_single_point_index()
        if idx is None:
            return

        pass

    def on_make_geodesic(self):
        pass

    def on_flip_out(self):
        pass

    def on_flip_edge(self):
        e = self.path_picker.get_corresponding_edge(self.tri)
        if e is None:
            return

        self.path_picker.on_clear()

        self.tri.flip(e)
        self.set_up_intrinsic_mesh()
        self.plotter.remove_actor('tri')
        tri_kwargs = {
            'name': 'tri',
            'edge_color': 'Black',
            'line_width': TRIANGULATION_EDGE_WIDTH,
            'show_edges': SHOW_TRIANGULATION_EDGES,
            'show_scalar_bar': False,
            'cmap': TRIANGULATION_FACES_COLOR_MAP
        }
        self.tri_obj.set_active_scalars('TriagColor')

        self.plotter.add_mesh(self.tri_obj, **tri_kwargs)




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

