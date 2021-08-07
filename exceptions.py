from math import degrees


class TriangulationException(BaseException):
    pass


class NonExistentEdgeException(TriangulationException):
    def __init__(self, origin, dst):
        msg = "Edge {}->{}, " \
              "does not exist" \
            .format(origin, dst)
        super().__init__(msg)


class NonExistentJointException(TriangulationException):
    def __init__(self, a, b, c):
        msg = "The vertices {}->{}->{} " \
              "do not form a path along the intrinsic triangulation edges" \
            .format(a, b, c)
        super().__init__(msg)


class MistriangulationException(TriangulationException):
    def __init__(self, e):
        msg = "Cannot flip edge {}->{}, " \
              "because on of its incident faces are not a triangle" \
            .format(e.origin, e.dst)
        super().__init__(msg)


class LowDegreeVertexException(TriangulationException):
    def __init__(self, e, idx, deg):
        msg = "Cannot flip edge {}->{}, " \
              "because vertex {} has degree {}" \
            .format(e.origin, e.dst,
                    idx, deg)
        super().__init__(msg)


class LowDegreeVertexWarning(TriangulationException):
    def __init__(self, e, idx, deg):
        msg = "Flipping edge {}->{}, " \
              "despite vertex {} having degree {}" \
            .format(e.origin, e.dst,
                    idx, deg)
        super().__init__(msg)


class ReflexAngleException(TriangulationException):
    def __init__(self, e):
        msg = "Cannot flip edge {}->{}, " \
              "due to a reflex angle" \
            .format(e.origin, e.dst)
        super().__init__(msg)


class WedgeReflexAngleException(TriangulationException):
    def __init__(self, angle):
        msg = "The wedge angle {:.2f}° is reflex"\
            .format(degrees(angle))
        super().__init__(msg)


class SelfEdgeException(TriangulationException):
    def __init__(self, a, b, c):
        msg = "Cannot flip-out joint {}->{}->{}, " \
              "because its edges are not distinct" \
            .format(a, b, c)
        super().__init__(msg)


class IntersectionNotFoundException(TriangulationException):
    def __init__(self):
        msg = "Intersection not found"
        super().__init__(msg)
