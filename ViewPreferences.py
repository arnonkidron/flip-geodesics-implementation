# edge flips
REFLEX_ANGLE_THRESHOLD = 1e-02

# mesh edges vs. triangulation edges
MESH_EDGE_WIDTH = 2
TRIANGULATION_EDGE_WIDTH = 4
MESH_EDGE_COLOR = 'Black'
TRIANGULATION_EDGE_COLOR = MESH_EDGE_COLOR
SHOW_MESH_EDGES = True
SHOW_TRIANGULATION_EDGES = True

# triangulation faces
SHOW_TRIANGULATION_FACES = True
# if True:
TRIANGULATION_FACES_COLOR_MAP = 'Accent'
TRIANGULATION_FACES_COLOR_MAP_SIZE = 8
# if False:
MESH_FACE_COLOR = 'Gold'

# picked path
PICKED_POINT_SIZE = 50
PICKED_PATH_WIDTH = PICKED_POINT_SIZE / 2
PICKED_PATH_COLOR = 'Pink'
SHOW_PICKED_PATH = True
SHOW_PICKED_PATH_END_POINTS = True
SHOW_PICKED_PATH_ALL_POINTS = False

# picking the intersection points of triangulation edges with mesh edges
ALLOW_PICK_INTERSECTION_POINTS = True
PICKED_INTERSECTION_POINTS_COLOR = 'Red'

# compute the intersection points one at a time
COMPUTE_INTERSECTION_POINTS_ONE_AT_A_TIME = True
COMPUTED_INTERSECTION_POINTS_COLOR = PICKED_INTERSECTION_POINTS_COLOR
COMPUTED_INTERSECTING_EDGE_COLOR = PICKED_INTERSECTION_POINTS_COLOR
KEY_EVENT_RE_RENDER = 'a'

# key bindings
KEY_EVENT_CLEAR_PICKED_PATH = 'c'
KEY_EVENT_UNDO_PICK = 'z'
KEY_EVENT_SHOW_INFO = 'l'
KEY_EVENT_PICK_NEXT_EDGE = 'n'
KEY_EVENT_PICK_TWIN_EDGE = 'm'
KEY_EVENT_EDGE_FLIP = 'u'
KEY_EVENT_FLIPOUT = 'i'
KEY_EVENT_MAKE_GEODESIC = 'o'
KEY_EVENT_SINGLE_SOURCE_DIJKSTRA = '['
