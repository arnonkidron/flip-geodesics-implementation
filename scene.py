import pyvista as pv
from tkinter import Tk, messagebox, simpledialog
from tkinter.filedialog import askopenfilename
from multiple_path_picker import *
from path_shortener import *
from single_src import *
import view_preferences as prefer


class Scene:
    def __init__(self, url=None):
        if url is None:
            url = self.ask_for_url()

        self.V = None

        self.plotter = pv.Plotter()

        self.shrank_mesh_face_actor = None
        self.mesh_actor = self.set_up_extrinsic_mesh(url)
        self.tri = self.set_up_triangulation()
        self.tri_edge_actor = None
        self.tri_face_actor = None

        self.slow_edge = None
        self.slow_generator = None
        self.slow_generator2 = None
        self.result_path = MultiplePathVisualizer(self)

        self.add_extrinsic_mesh_actor()
        self.add_intrinsic_triangulation()

        self.path_picker = MultiplePathPicker(self)
        self.text_actor = None
        self.set_up_events()

    @staticmethod
    def ask_for_url():
        # we don't want a full GUI, so keep the root window from appearing
        Tk().withdraw()

        # show an "Open" dialog box and return the path to the selected file
        mesh_url = askopenfilename(
            title='Choose a mesh file to load',
            filetypes=[
                ('All Files (*.*)', '*'),
            ]
        )

        if mesh_url == '':
            exit("flip_geodesic load error: No file chosen")

        return mesh_url

    @staticmethod
    def ask_for_integer(msg):
        tk = Tk()
        tk.withdraw()
        answer = simpledialog.askstring(title="Input", prompt=msg, parent=tk)
        if answer is None:
            return None
        return int(answer)

    @staticmethod
    def warn(title, msg):
        """
        Outcome: display the message in a tkinter window, with retry & cancel
        Output: whether the user has chosen to retry
        """
        print(msg)
        Tk().withdraw()
        return messagebox.showwarning(title=title, message=msg, type=messagebox.OK)

    def show(self):
        self.plotter.show()

    def set_up_extrinsic_mesh(self, url):
        obj = pv.read(url)

        if hasattr(obj, 'is_all_triangles'):
            if not obj.is_all_triangles():
                obj.triangulate()

        V = obj.points
        if hasattr(obj, 'faces'):
            E = obj.faces
        else:
            E = obj.cells

        self.mesh_actor = pv.PolyData(V, lines=E)

        obj.compute_normals(inplace=True)

        self.shrank_mesh_face_actor = obj.warp_by_vector(factor=prefer.UNDERLYING_SHRINK_FACTOR)

        return self.mesh_actor

    def set_up_triangulation(self):
        V = self.mesh_actor.points
        F = self.mesh_actor.lines
        F = np.reshape(F, (len(F) // 4, 4))
        F = np.delete(F, 0, axis=1)
        tri = IntrinsicTriangulation(V, F)
        tri.init_coloring(len(F), prefer.TRIANGULATION_FACES_COLOR_MAP_SIZE)

        self.tri = tri
        return tri

    def add_extrinsic_mesh_actor(self):
        mesh_kwargs = {
            'name': 'extrinsic_mesh',
            'edge_color': prefer.MESH_EDGE_COLOR,
            'color': prefer.MESH_EDGE_COLOR,
            'line_width': prefer.MESH_EDGE_WIDTH,
            'show_edges': prefer.SHOW_MESH_EDGES,
            'show_scalar_bar': False,
        }

        if prefer.SHOW_TRIANGULATION_FACES:
            mesh_kwargs['style'] = 'wireframe'
        else:
            mesh_kwargs['color'] = prefer.MESH_FACE_COLOR

        if prefer.LABEL_VERTICES:
            self.plotter.add_point_labels(self.mesh_actor, range(self.mesh_actor.n_points),
                                           font_size=10)

        # self.plotter.remove_actor(mesh_kwargs['name'])
        self.plotter.add_mesh(self.mesh_actor, **mesh_kwargs)

        if prefer.SHOW_UNDERLYING_MESH_FACES:
            shrank_mesh_kwargs = {
                'name': 'shrank_extrinsic_mesh',
                'color': prefer.UNDERLYING_COLOR,
                'show_edges': False,
                'show_scalar_bar': False,
            }
            self.plotter.add_mesh(self.shrank_mesh_face_actor, **shrank_mesh_kwargs)

    def add_intrinsic_triangulation(self):
        # advance intersection point computation by one step
        if self.slow_generator is not None:

            intersection = next(self.slow_generator)
            if intersection is None:
                self.slow_edge = None
                self.slow_generator = None
                self.plotter.remove_actor('slow_edge')
                self.plotter.remove_actor('hit_edge')
                self.plotter.remove_actor('hit_point')
                self.plotter.remove_actor('next_vec')
            else:
                point, mesh_edge = intersection.coords, intersection.mesh_edge
                actor = pv.PolyData([self.mesh_actor.points[self.slow_edge.origin], self.mesh_actor.points[self.slow_edge.dst]], [2, 0, 1])
                self.plotter.add_mesh(actor, name='slow_edge', style='wireframe', color=prefer.RESULT_PATH_COLOR, render_lines_as_tubes=True, line_width=prefer.PATH_EDGE_WIDTH)
                actor = pv.PolyData([self.mesh_actor.points[mesh_edge.origin], self.mesh_actor.points[mesh_edge.dst]], [2, 0, 1])
                self.plotter.add_mesh(actor, name='hit_edge', style='wireframe', color=prefer.COMPUTED_INTERSECTING_EDGE_COLOR, render_lines_as_tubes=True, line_width=prefer.PATH_EDGE_WIDTH)
                self.plotter.add_mesh(pv.PolyData(point), name='hit_point', color=prefer.COMPUTED_INTERSECTION_POINTS_COLOR, render_points_as_spheres=True, point_size=prefer.PATH_POINT_SIZE)
                actor = pv.PolyData()
                for out_vec in intersection.get_out_vectors():
                    actor += pv.Arrow(intersection.coords, out_vec, tip_radius=0.25, shaft_radius=0.10, scale='auto')
                self.plotter.add_mesh(actor, name='next_vec', color='Green')

        # set up
        self.V, E, F, coloring = self.tri.get_poly_data(self.tri.mesh, need_extrinsic_faces=prefer.SHOW_TRIANGULATION_FACES)

        self.tri_edge_actor = pv.PolyData(self.V, E)
        self.tri_face_actor = pv.PolyData(self.V, F)
        self.tri_face_actor.cell_arrays['TriagColoring'] = coloring
        self.tri_face_actor.set_active_scalars('TriagColoring')

        tri_edge_kwargs = {
            'name': 'tri_edge',
            'style': 'wireframe',
            'color': prefer.TRIANGULATION_EDGE_COLOR,
            'line_width': prefer.TRIANGULATION_EDGE_WIDTH,
            'show_scalar_bar': False,
        }

        tri_face_kwargs = {
            'name': 'tri_face',
            'show_edges': False,
            'show_scalar_bar': False,
            'cmap': prefer.TRIANGULATION_FACES_COLOR_MAP,
        }

        # self.plotter.remove_actor(tri_edge_kwargs['name'])
        # self.plotter.remove_actor(tri_face_kwargs['name'])
        if prefer.SHOW_TRIANGULATION_EDGES:
            self.plotter.add_mesh(self.tri_edge_actor, **tri_edge_kwargs)

        if prefer.SHOW_TRIANGULATION_FACES:
            self.plotter.add_mesh(self.tri_face_actor, **tri_face_kwargs)

        self.result_path.add_actor()

    def set_up_events(self):
        bindings = {
            prefer.KEY_EVENT_MAKE_GEODESIC: self.on_make_geodesic,
            prefer.KEY_EVENT_FLIPOUT: self.on_flip_out,
            prefer.KEY_EVENT_EDGE_FLIP: self.on_edge_flip,
            prefer.KEY_EVENT_DELAUNAY: self.on_delaunay,
            prefer.KEY_EVENT_SHOW_INFO: self.on_info,
            prefer.KEY_EVENT_CLEAR_PICKED_PATH: self.on_clear,
            prefer.KEY_EVENT_START_NEW_PATH: self.path_picker.on_start_new_path,
            prefer.KEY_EVENT_UNDO_PICK: self.path_picker.on_undo,
            prefer.KEY_EVENT_CLOSE_LOOP: self.path_picker.on_close_loop,
            prefer.KEY_EVENT_PICK_NEXT_EDGE: self.path_picker.on_pick_next_edge,
            prefer.KEY_EVENT_PICK_TWIN_EDGE: self.path_picker.on_pick_twin_edge,
            prefer.KEY_EVENT_RE_RENDER: self.add_intrinsic_triangulation,
            prefer.KEY_EVENT_PICK_RESULT_PATH: self.on_pick_result_path,
            prefer.KEY_EVENT_CLEAR_RESULT_PATH: self.on_clear_result_path,
            prefer.KEY_EVENT_PICK_BY_INDEX: self.on_pick_by_index,
            prefer.KEY_EVENT_SHOW_VECS_ONE_AT_A_TIME: self.on_show_vecs_one_at_a_time,
        }
        for key, callback in bindings.items():
            self.plotter.add_key_event(key, callback)

    def on_pick_by_index(self, idx=None):
        if idx is None:
            idx = self.ask_for_integer("Enter a vertex index to pick:")
            if idx is None:
                return

        self.path_picker.on_pick(self.V[idx])

    def on_single_source(self):
        idx = self.path_picker.get_single_point_index()
        if idx is None:
            return

        kwargs = {
            # 'color': prefer.TRIANGULATION_EDGE_COLOR,
            'line_width': prefer.TRIANGULATION_EDGE_WIDTH,
            'point_size': 0,
        }

        shortener = SingleSrcShortener(self.tri)
        path = self.path_picker.get_path()
        shortener.set_path(deepcopy(path))
        try:
            shortener.make_geodesic()
            print("Finished computing single source geodesics successfully. Start rendering")
        except TriangulationException as err:
            self.result_path.set_path(shortener.get_path(), **kwargs)
            self.result_path.set_path(shortener.get_frontier(), **kwargs, color='Black')
            self.add_intrinsic_triangulation()
            return self.warn(title="MakeGeodesic fail", msg=str(err))

        self.result_path.set_path(shortener.get_path(), **kwargs)
        self.add_intrinsic_triangulation()

    def on_delaunay(self):
        excluded_edges = self.path_picker.get_path_edge_tuples_set().union(
            self.result_path.get_path_edge_tuples_set()
        )

        self.tri.delaunay(excluded_edges)
        self.add_intrinsic_triangulation()

    def on_make_geodesic(self):
        roi = self.path_picker.roi
        if roi == ROI.VERTEX:
            return self.on_single_source()
        elif roi == ROI.EMPTY or roi == ROI.INTERSECTION:
            return

        shortener = get_shortener(roi, self.tri)
        path = self.path_picker.get_path()

        try:
            shortener.set_path(deepcopy(path))
        except NonExistentJointException:
            # cannot redo
            return
        try:
            shortener.make_geodesic()
            new_path = shortener.get_path()
            self.set_result(new_path, prefer.SHOW_ON_MAKE_GEODESIC)
        except TriangulationException as err:
            self.add_intrinsic_triangulation()
            return self.warn(title="MakeGeodesic fail", msg=str(err))

    def on_flip_out(self):
        if self.path_picker.roi == ROI.VERTEX:
            if self.slow_generator2 is None:
                shortener = SingleSrcShortener(self.tri)
                shortener.set_path(self.path_picker.get_path())
                self.slow_generator2 = shortener.make_geodesic_one_at_a_time()
            next(self.slow_generator2)
            self.add_intrinsic_triangulation()
            return

        if self.result_path.is_empty():
            viz = self.path_picker
        else:
            viz = self.result_path
        shortener = get_shortener(viz.roi, self.tri)
        if shortener is None:
            return

        shortener.set_path(deepcopy(viz.get_path()))
        try:
            shortener.flipout_the_minimal_wedge()
            new_path = shortener.get_path()
            self.set_result(new_path, prefer.SHOW_ON_FLIPOUT)
        except TriangulationException as err:
            self.add_intrinsic_triangulation()
            return self.warn(title="FlipOut fail", msg=str(err))

    def on_edge_flip(self):
        if self.path_picker.roi != ROI.EDGE:
            return

        old_edge = self.path_picker.get_corresponding_edge()
        try:
            new_edge = self.tri.flip(old_edge)
        except TriangulationException as err:
            self.add_intrinsic_triangulation()
            return self.warn(title="Edge flip fail", msg=str(err))

        # compute intersection points
        if prefer.COMPUTE_INTERSECTION_POINTS_ONE_AT_A_TIME:
            self.slow_edge = new_edge
            self.slow_generator = new_edge.init_intersections_one_at_a_time(self.tri.mesh)
        else:
            new_edge.init_intersections(self.tri.mesh)

        self.set_result([new_edge.origin, new_edge.dst], prefer.SHOW_ON_EDGE_FLIP)
        self.remove_text()

    def set_result(self, path, what_to_show):
        if what_to_show == prefer.Show.NOTHING:
            return self.on_clear()

        self.result_path.set_path(path)

        if what_to_show == prefer.Show.ONLY_THE_RESULT:
            self.on_pick_result_path()

        self.add_intrinsic_triangulation()

    def on_info(self):
        self.remove_text()

        msg = None
        roi = self.path_picker.roi
        if roi == ROI.VERTEX:
            idx = self.path_picker.get_single_point_index()
            coords = self.V[idx]
            deg = len(self.tri.in_edges[idx])
            msg = "Vertex #{}\n" \
                  "({:.4f},{:.4f},{:.4f})\n" \
                  "Degree {}\n"\
                .format(idx, coords[0], coords[1], coords[2], deg)
        elif roi == ROI.INTERSECTION:
            idx = self.path_picker.get_single_point_index()
            coords = self.V[idx]
            msg = "Intersection point\n" \
                  "({:.4f},{:.4f},{:.4f})\n"\
                .format(coords[0], coords[1], coords[2])
        elif roi == ROI.EDGE:
            e = self.path_picker.get_corresponding_edge()
            msg = e.get_info()
            self.show_edge_first_vec(e)
        else:
            msg = ""
            for path in self.path_picker.get_paths_for_info():
                if len(path) == 0:
                    continue
                if path[0] == path[-1]:
                    msg += "Loop "
                else:
                    msg += "Path "
                for v in path:
                    msg += str(v)
                    msg += "->"

                msg = msg[:-2]
                msg += "\n"

        if msg is None:
            self.text_actor = None
        else:
            self.text_actor = self.plotter.add_text(msg, name='info')
            print(msg)

    def on_show_vecs(self):
        idx = self.path_picker.get_single_point_index()
        if idx is not None:
            self.show_vec_from_vertex(idx)

        e = self.path_picker.get_corresponding_edge()
        if e is not None:
            self.show_vec_from_edge(e)

    def show_vec_from_vertex(self, idx):
        coords = self.V[idx]
        actor = pv.PolyData()
        for in_edge in self.tri.in_edges[idx]:
            out_edge = in_edge.twin
            vec = out_edge.get_first_segment_vector()

            actor += pv.Arrow(coords, vec, tip_radius=0.25, shaft_radius=0.10, scale='auto')

        self.plotter.add_mesh(actor, name='vecs', color='SpringGreen')

    def show_edge_first_vec(self, e):
        actor = pv.PolyData()
        actor += pv.Arrow(self.V[e.origin], e.get_first_segment_vector(), tip_radius=0.25, shaft_radius=0.10, scale='auto')
        self.plotter.add_mesh(actor, name='vecs', color='SpringGreen')

    def show_vec_from_edge(self, e):
        actor = pv.PolyData()
        actor += pv.Arrow(self.V[e.origin], e.get_first_segment_vector(), tip_radius=0.25, shaft_radius=0.10, scale='auto')

        e.init_intersections(self.tri.mesh)
        if e.intersections_status != e.Status.FAILED:
            for intersection in e.get_intersections(self.tri.mesh):
                for out_vec in intersection.get_out_vectors():
                    actor += pv.Arrow(intersection.coords, out_vec, tip_radius=0.25, shaft_radius=0.10, scale='auto')
        else:
            pass

        self.plotter.add_mesh(actor, name='vecs', color='SpringGreen')

    def on_show_vecs_one_at_a_time(self):
        idx = self.path_picker.get_single_point_index()
        if idx is not None:
            self.show_vec_from_vertex(idx)

        if self.slow_edge is None:
            self.slow_generator = self.show_vecs_one_at_a_time()
            self.slow_edge = next(self.slow_generator)
            if self.slow_edge is None:
                self.warn("", "No single edge selected")
                return

        next(self.slow_generator)

    def show_vecs_one_at_a_time(self):
        e = self.path_picker.get_corresponding_edge()
        yield e
        if e is None:
            return

        arrows = pv.Arrow(self.V[e.origin], e.intersections_.get_first_segment_vector(), tip_radius=0.25, shaft_radius=0.10, scale='auto')
        # self.plotter.remove_actor('vecs')
        self.plotter.add_mesh(arrows, name='vecs', color='SpringGreen')
        yield True

        for intersection in e.get_intersections(self.tri.mesh):
            if intersection.is_fake:
                break

            point = pv.PolyData(intersection.coords)
            self.plotter.add_mesh(point, name='intersection', color=prefer.PICKED_INTERSECTION_POINTS_COLOR, render_points_as_spheres=True, point_size=prefer.PATH_POINT_SIZE)

            arrows = pv.PolyData()
            for out_vec in intersection.get_out_vectors():
                arrows += pv.Arrow(intersection.coords, out_vec, tip_radius=0.25, shaft_radius=0.10, scale='auto')
            self.plotter.add_mesh(arrows, name='vecs', color='SpringGreen')
            yield True

        self.plotter.remove_actor('intersection')
        self.plotter.remove_actor('vecs')
        self.slow_generator = None
        self.slow_edge = None
        yield False

    def remove_text(self):
        self.plotter.remove_actor(self.text_actor)
        self.plotter.remove_actor('vecs')
        self.slow_generator = None
        self.slow_edge = None

    def on_clear(self):
        self.path_picker.clear()
        self.result_path.clear()
        self.remove_text()

    def on_pick_result_path(self):
        if self.result_path.is_empty():
            return

        self.path_picker.set_path(self.result_path.get_path())
        self.result_path.clear()

    def on_clear_result_path(self):
        self.result_path.clear()
