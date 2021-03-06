from math import radians

# thresholds
FLAT_ANGLE_THRESHOLD_FOR_EDGE_FLIP = radians(1)
FLAT_ANGLE_THRESHOLD_FOR_FLIP_OUT = radians(1)
REFLEX_ANGLE_THRESHOLD_FOR_SINGLE_SRC_GEODESICS = 1e-09

# limits for make_geodesic
ITERATIONS_LIMIT = None
LENGTH_THRESHOLD = None
