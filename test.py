import polyscope as ps
from main import init_polyscope
from utils import *
from math import radians
import numpy as np


def test_rotate():
    origin = (0,0,0)
    a = (1,0,0)
    b = (0,1,0)
    c = turn(a, radians(22.5), towards=b)
    d = turn(c, radians(22.5), towards=b)
    e = turn(b, radians(22.5), towards=d)

    V = np.array([origin, a, b, c, d, e])

    ps = init_polyscope()
    ps.register_curve_network("a", V, np.array([[0,1]]))
    ps.register_curve_network("b", V, np.array([[0,2]]))
    ps.register_curve_network("c", V, np.array([[0,3]]))
    ps.register_curve_network("d", V, np.array([[0,4]]))
    ps.register_curve_network("e", V, np.array([[0,5]]))
    ps.show()


test_rotate()
