from edge import *
from Triangulation import *
import numpy as np
from utils import is_reflex

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
        wedge_angle = self.tri.get_wedge_angle(a, b, c)
        if is_reflex(wedge_angle):
            print("The wedge angle {:.2f}° is not reflex".format(degrees(wedge_angle)))
            return

        # find edges of this wedge
        e1 = self.tri.get_edge(a, b)
        e2 = self.tri.get_edge(b, c)
        if e1 is None or e2 is None:
            return

        # calculate bypassing path
        bypass = [a]
        e = e1.next
        while e.dst != c:
            bypass.append(e.dst)
            e = e.twin.next
        bypass.append(c)

        i = 1
        while i != len(bypass) - 1:
            pass
            # if fail:
            #     i = i - 1
            # else:
            #     i = i + 1

        return bypass

        # update self.path
        # remove b
        # add bypass instead of it

    def flipout_the_minimal_wedge(self):
        # find minimal wedge
        if len(self.path) != 3:
            return "We'll program this later"
        a, b, c = self.path[0], self.path[1], self.path[2]

        # call self.flipout()
        return self.flipout(a, b, c)

    def make_geodesic(self):
        while not self.is_geodesic() and False:
            self.flipout_the_minimal_wedge()

    def is_geodesic(self):
        return False
