from exceptions import NonExistentJointException
from enum import Enum


class Joint:
    """
    Joints are paths that consist of 2 edges, a->b->c.
    They are the building blocks of paths for FlipOut.

    :field e1, e2: the edges composing the joint. e1 is a->b, and e2 is b->c
    :field path_data: to be used by PathShortener to store anything it needs,
                      e.g. the index of vertex b in the path
    :field wedge_angle: the angle to the left of the joint
    :field direction: if the joint belongs to a path ...->a->b->c->..., then it goes forth
                     otherwise, it belongs to a path ...->c->b->a->..., so it goes back
    :field twin: the joint c->b->a
    """
    def __init__(self, edge_1, edge_2, path_data, direction, twin=None):
        self.e1 = edge_1
        self.e2 = edge_2
        self.length = edge_1.length + edge_2.length
        self.wedge_angle = self.init_wedge_angle()

        self.path_data = path_data

        self.direction = direction
        self.twin = twin

    class Direction(Enum):
        FORTH = 1
        BACK = -1

    def init_wedge_angle(self):
        if self.e1 == self.e2:
            return 0

        e = self.e1.next
        sum = 0
        while e != self.e2:
            sum += e.corner_angle
            e = e.twin.next
        sum += e.corner_angle

        self.wedge_angle = sum
        return sum

    @staticmethod
    def from_vertices(triangulation, a, b, c, path_data):
        edge_1 = triangulation.get_edge(a, b)
        edge_2 = triangulation.get_edge(b, c)
        
        if edge_1 is None or edge_2 is None:
            raise NonExistentJointException(a, b, c)
        
        joint = Joint(edge_1, edge_2, path_data, Joint.Direction.FORTH)
        twin = Joint(edge_2.twin, edge_1.twin, path_data, Joint.Direction.BACK, twin=joint)
        return joint

    @property
    def twin(self):
        return self._twin

    @twin.setter
    def twin(self, twin):
        self._twin = twin
        if twin is not None:
            twin._twin = self

    @property
    def vertices(self):
        return \
            self.e1.origin, \
            self.e2.origin, \
            self.e2.dst
