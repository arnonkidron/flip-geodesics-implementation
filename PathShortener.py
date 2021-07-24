from edge import *
from Triangulation import *
import numpy as np


class PathShortener:
    def __init__(self, triangulation):
        self.tri = triangulation
        self.path = []
        # [543, 65, 47, 82]
        # [543, 65, 82]
        # [543, 65, 64, 24, 35, 82]

    def set_path(self, path):
        self.path = path
        if len(path) > 2:
            wedge_angles = [self.tri.get_wedge_angle(
                self.path[i],
                self.path[i+1],
                self.path[i+2],
            ) for i in range(len(path) - 2)]

    def get_path(self):
        return self.path

    def flipout(self, a, b, c):
        pass

    def flipout_the_minimal_wedge(self):
        print("Hi")
        # find minimal wedge
        # call self.flipout()
        pass

    def make_geodesic(self):
        while not self.is_geodesic() and False:
            self.flipout_the_minimal_wedge()

    def is_geodesic(self):
        return False
