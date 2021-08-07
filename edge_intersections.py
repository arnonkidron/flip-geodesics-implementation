import numpy as np
from utils import rotate, turn
from math import pi
from enum import Enum
import view_preferences as prefer
from copy import copy


class Intersection:
    def __init__(self, coords, mesh_edge=None, distance_from_mesh_edge=None, in_vec=None, out_vec=None, is_fake=False):
        self.coords = coords
        self.mesh_edge = mesh_edge
        self.distance_from_mesh_edge = distance_from_mesh_edge
        self.distance_from_face_center = None if self.mesh_edge is None else \
            np.linalg.norm(self.coords - self.mesh_edge.get_face_normal())
        self.in_vec = in_vec
        self.out_vec = out_vec
        self.is_fake = is_fake
        self.other_candidates = []

    @property
    def error(self):
        return self.distance_from_face_center

    @property
    def surface_error(self):
        n = self.mesh_edge.get_face_normal()
        return np.dot(n, self.in_vec)

    def get_next_edges(self):
        return [
            self.mesh_edge.twin.next,
            self.mesh_edge.twin.next.next,
        ]

    def get_out_vectors(self):
        if self.out_vec is not None:
            return [self.out_vec]

        if self.mesh_edge is None or self.in_vec is None:
            return []

        out_vec1 = rotate(self.in_vec, self.mesh_edge.mesh_face_angle, self.mesh_edge.vec)
        out_vec2 = rotate(self.in_vec, 2 * pi - self.mesh_edge.mesh_face_angle, self.mesh_edge.vec)
        return [out_vec1, out_vec2]

    def test_is_vec_on_face(self):
        n = self.mesh_edge.get_face_normal()
        dot_products = [np.dot(vec / np.linalg.norm(vec), n) for vec in self.get_out_vectors()]
        print(dot_products)


class EdgeIntersections:
    def __init__(self, edge):
        self.edge = edge
        self.status = self.Status.UNINITIALIZED
        self.points = []

        self.best_attempt_points = None
        self.best_attempt_score = np.inf

    class Status(Enum):
        UNINITIALIZED = 1
        FINISHED = 2
        TWIN_FINISHED = 3
        FAILED = 4

    @property
    def twin(self):
        return self.edge.twin.intersections_

    def set_no_intersection(self):
        self.status = self.Status.FINISHED
        self.points = []

    def get_start_point(self):
        return Intersection(self.edge.near_mesh_edge.origin_coords)

    def get_first_segment_vector(self):
        next_mesh_edge = self.edge.near_mesh_edge.twin.next
        return turn(
            self.edge.near_mesh_edge.vec,
            self.edge.angle_from_near_mesh_edge,
            towards=next_mesh_edge.vec
        )

    def get_first_intersecting_mesh_edge(self):
        return self.edge.near_mesh_edge.twin.next.next

    def compute(self, mesh):
        for _ in self.compute_one_at_a_time(mesh):
            pass

    def compute_one_at_a_time(self, mesh):
        if self.status != self.Status.UNINITIALIZED:
            yield None; return

        self.points = []
        dst = self.edge.dst

        if self.edge.near_mesh_edge.is_point_in_face(dst):
            print("--------Superfluous---------")
            self.edge.sit_on_mesh_edge(self.edge.near_mesh_edge)
            yield None; return

        # initial values
        prev_intersection = self.get_start_point()
        vecs = [self.get_first_segment_vector()]
        mesh_edges = [self.get_first_intersecting_mesh_edge()]

        while len(self.points) < 200:
            # try both edges, and both vectors
            intersection = None
            candidates = []
            for e in mesh_edges:
                for vec in vecs:
                    candidate = e.get_intersection(prev_intersection.coords, vec)
                    if candidate is not None:
                        candidates.append(candidate)

            if len(candidates) == 1:
                intersection = candidates[0]
            elif len(candidates) > 1:
                candidates.sort(key=lambda x: x.distance_from_mesh_edge)
                intersection = candidates[0]
                intersection.other_candidates = candidates[1:]
            elif len(candidates) == 0:
                intersection = None
                # intersection = mesh.get_intersection_complete_search(prev_intersection.coords, vecs, mesh_edges[0])
                # if intersection is not None \
                #         and intersection.distance_from_mesh_edge > INTERSECTION_THRESHOLD:
                #     intersection = None

            if intersection is None:
                self.status = self.Status.FAILED

                # let the twin try his luck
                self.twin.compute(mesh)
                if self.twin.status == self.Status.FINISHED:
                    # if he succeeds - no failure, take his intersections
                    self.status = self.Status.TWIN_FINISHED
                    self.points = []
                    yield None; return
                else:
                    prev_intersection, intersection = self.backtrack(mesh)
                    if prev_intersection is not None:
                        continue
                    yield None; return

            self.points.append(intersection)

            vecs = intersection.get_out_vectors()
            mesh_edges = intersection.get_next_edges()

            prev_intersection.out_vec = intersection.in_vec
            prev_intersection = intersection

            yield intersection
            if mesh_edges[0].is_point_in_face(dst):
                break

        self.status = self.Status.FINISHED
        if self.twin.status == self.Status.UNINITIALIZED:
            self.twin.status = self.Status.TWIN_FINISHED
        yield None; return

    def get_intersections(self, mesh):
        if self.status == self.Status.UNINITIALIZED:
            self.compute(mesh)

        if self.status == self.Status.TWIN_FINISHED:
            return np.flipud(self.twin.get_intersections(mesh))

        return self.points

    def backtrack(self, mesh):
        # we couldn't compute the intersection points properly. Now we go back and try other candidates

        # save attempt
        # score = sum([intersection.surface_error for intersection in self.points])
        score = sum([np.linalg.norm(self.points[i].coords - self.points[i+1].coords) for i in range(len(self.points) - 1)])
        if self.best_attempt_score > score:
            self.best_attempt_score = score
            self.best_attempt_points = copy(self.points)

        # backtrack
        while len(self.points) > 1:
            self.points.pop()
            other_candidates = self.points[-1].other_candidates
            if len(other_candidates) > 0:
                self.points.pop()
                intersection = other_candidates[0]
                intersection.other_candidates = other_candidates[1:]
                prev_intersection = self.points[-1] if len(self.points) > 0 else self.get_start_point()
                return prev_intersection, intersection

        # end of backtracking
        self.edge.print("---------- edge", "failed to compute intersection points-----------")
        self.points = self.best_attempt_points
        self.wrap_up_failure(mesh)
        return None, None

    def wrap_up_failure(self, mesh):
        # we failed to compute all intersection points, and now we need to make it look nice for rendering
        # add another point, lifted above the previous one, in order to make the path go above the mesh
        last_intersection = self.points[-1]
        e = last_intersection.mesh_edge
        normal = e.get_face_normal() + e.twin.get_face_normal()
        normal /= np.linalg.norm(normal)
        normal *= mesh.avg_face_area
        normal *= prefer.FAILED_EDGES_JUMP_HEIGHT_COEF
        lift_coords = last_intersection.coords + normal
        lift = Intersection(lift_coords, is_fake=True)
        self.points.append(lift)
