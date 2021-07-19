import polyscope as ps
import numpy as np
from FlipEdgeNetwork import FlipEdgeNetwork

def init_polyscope():
    # enable auto centering and scaling
    # ps.set_autocenter_structures(True)
    # ps.set_autoscale_structures(True)

    # Camera ctrls
    ps.set_navigation_style("free")
    ps.set_up_dir("z_up")

    # initialize
    ps.set_program_name("Flip Geodesics Implementation for 236329 GDP")
    ps.set_verbosity(2)
    ps.set_use_prefs_file(False)
    ps.set_errors_throw_exceptions(True)
    ps.init()

    return ps

from tkinter import Tk     # from tkinter import Tk for Python 3.x
from tkinter.filedialog import askopenfilename
def open_mesh_file():
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    mesh_url = askopenfilename() # show an "Open" dialog box and return the path to the selected file
    if mesh_url == '':
        exit("No file chosen")
    return open(mesh_url, "r")

def get_example_shape():
    V = np.array([ [0, 5., 0], [0, 1, -3.], [-4., 0, 0], [0, 1, 3.], [4., 0, 0] ])
    F = np.array([ [0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1], [1, 4, 2], [2, 4, 3] ])
    return V, F


# mesh_file = open_mesh_file()
# V,F = parse_mesh_file(mesh_file)

V, F = get_example_shape()
net = FlipEdgeNetwork(V, F)
net.set_path([0,1,2])
pathV, pathE = net.get_path_polyline()
triV, triE = net.get_intrinsic_trinagulation_polyline()

ps = init_polyscope()

ps_mesh = ps.register_surface_mesh("Mesh", V, F)
ps_path = ps.register_curve_network("Path", pathV, pathE)
# ps_tri = ps.register_curve_network("Intrinsic Triangulation", triV, triE)


# ps_e = ps.register_curve_network("Path", V, np.array([[0, 1], [1,2]]))
# ps_t = ps.register_curve_network("Triangulation Edges", V, np.array([[0, 1], [1,2]]))
# ps_v = ps.register_point_cloud("Path Vertices", np.array([V[0],V[1],V[2]]))

## generate some random nodes and edges between them
# nodes = np.random.rand(100, 3)
# edges = np.random.randint(0, 100, size=(250,2))
# ps_net = ps.register_curve_network("Mesh Edges", nodes, edges)


# Pass control flow to polyscope, displaying the interactive window.
# Function will return when user closes the window.
ps.show()


# ps.reset_camera_to_home_view()
