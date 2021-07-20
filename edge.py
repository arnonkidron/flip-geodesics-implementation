import numpy as np


class Edge:
    def __init__(self, endpoint1, endpoint2):
        """

        :arg endpoint1,2: the indices of the two vertices connected by this edge.
        We require them in the ctor for the sake of future needs.
        Currently we do not save them because we rely on the matrix indices.

        :arg midpoints - the coordinates of the points that connect the
        different straight segments of this edge
        Note: the midpoint list is empty iff u,v have a mesh edge between them
        """
        self.endpoint1 = endpoint1
        self.endpoint2 = endpoint2
        self.midpoints = None
        self.length = None  #total length
        self.angle = None  # the angle between u &

        # line equation
        self.slope = None
        self.intercept = None

        self.midpoints = []  # TODO: remove

    def set_midpoints(self, lst):
        self.midpoints = lst

    def calc_midpoints(self, F, first_face_index, len2, len3):
        pass

    def get_polyline(self, index_difference=0):
        if self.midpoints is None:
            # possible error
            pass

        midpoint_indices = list(range(
            index_difference,
            len(self.midpoints) + index_difference
        ))

        verts = [self.endpoint1] \
            + midpoint_indices \
            + [self.endpoint2]

        return self.midpoints, [[verts[i], verts[i+1]] for i in range(len(verts) - 1)]

