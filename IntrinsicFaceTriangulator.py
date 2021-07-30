import numpy as np
import ViewPreferences as prefer


class IntrinsicFaceTriangulator:
    def __init__(self, e, intrinsic_face, mesh):
        self.intrinsic_face = intrinsic_face
        self.n = len(intrinsic_face)
        self.matched_vertices = []
        self.extrinsic_faces = []

        self.edges = [[SimpleHalfEdge(i - 1)] for i in range(self.n)]
        for i in range(self.n):
            self.edges[i-1][0].next = self.edges[i][0]

        uv = e
        vw = e.next
        wu = e.next.next
        u, v, w = uv.origin, vw.origin, wu.origin

        index_u = np.argmax(self.intrinsic_face == u)
        index_v = np.argmax(self.intrinsic_face == v)
        index_w = np.argmax(self.intrinsic_face == w)

        to_match_uv = self.put_intersections_to_match(uv.get_intersections(mesh), index_u, index_v, index_w)
        to_match_vw = self.put_intersections_to_match(vw.get_intersections(mesh), index_v, index_w, index_u)
        to_match_wu = self.put_intersections_to_match(wu.get_intersections(mesh), index_w, index_u, index_v)

        self.match_vertices(to_match_uv, [to_match_vw, to_match_wu])
        self.match_vertices(to_match_vw, [to_match_wu])

        self.compute_edges()
        self.compute_faces()

    def put_intersections_to_match(self, intersections, start, end, opposite):
        list_to_match = []
        if end == 0:
            end = self.n
        opposite_idx = self.intrinsic_face[opposite]
        for x in range(start + 1, end):
            e = intersections[x - start - 1].mesh_edge
            if opposite_idx in [e.origin, e.dst]:
                self.matched_vertices.append((x, opposite))
            else:
                list_to_match.append((x, e))

        return list_to_match

    def match_vertices(self, vertices_on_one_edge, other_edges_lists):
        for (x, e) in vertices_on_one_edge:
            for other_list in other_edges_lists:
                must_find = other_list == other_edges_lists[-1]
                y = self.find_and_delete_match(e, other_list, must_find)

                if y is not None:
                    self.matched_vertices.append((x, y))
                    break

    def find_and_delete_match(self, e, other_list, raise_not_found_exception=True):
        for i in range(len(other_list)):
            (y, some_e) = other_list[i]
            if some_e == e or some_e.twin == e:
                del other_list[i]
                return y

        if raise_not_found_exception \
                and not prefer.COMPUTE_INTERSECTION_POINTS_ONE_AT_A_TIME:
            raise "Rendering failure"

    def compute_edges(self):
        for (x, y) in self.matched_vertices:
            prev1 = self.edges[x][0]
            prev2 = self.edges[y][0]
            next1 = prev2.next
            next2 = prev1.next

            e1 = SimpleHalfEdge(x)
            e2 = SimpleHalfEdge(y)
            self.edges[y].append(e1)
            self.edges[x].append(e2)

            prev1.next = e1
            prev2.next = e2
            e1.next = next1
            e2.next = next2

    def compute_faces(self):
        for edge_list in self.edges:
            for e in edge_list:
                if e.visited:
                    continue

                face = [None]
                while not e.visited:
                    e.visited = True
                    vertex = self.intrinsic_face[e.origin]
                    face.append(vertex)
                    e = e.next

                face[0] = len(face) - 1
                self.extrinsic_faces.append(face)

    def get_faces(self):
        return np.hstack(self.extrinsic_faces), len(self.extrinsic_faces)


class SimpleHalfEdge:
    def __init__(self, origin):
        self.origin = origin
        self.next = None
        self.visited = False

