from PathPicker import *


class MultiplePathVisualizer:
    def __init__(self, scene):
        self.scene = scene
        self.visualizers = []

    @property
    def roi(self):
        if self.is_one():
            return self.last.roi
        else:
            return ROI.NETWORK

    def add_one(self, **kwargs):
        new_one = PathVisualizer(self.scene, name="_result_path_" + str(len(self.visualizers)), **kwargs)
        self.visualizers.append(new_one)

    def is_one(self):
        return len(self.visualizers) == 1

    def is_empty(self):
        return len(self.visualizers) == 0

    @property
    def last(self):
        return self.visualizers[-1]

    def get_path(self):
        if self.is_one():
            return self.last.get_path()
        else:
            return [x.get_path() for x in self.visualizers]

    def get_paths_for_info(self):
        return [x.indices for x in self.visualizers]

    def set_path(self, paths, **kwargs):
        if not isinstance(paths[0], list):
            paths = [paths]

        self.visualizers = []
        for path in paths:
            self.add_one(**kwargs)
            self.last.set_path(path)

    def add_actor(self):
        for x in self.visualizers:
            x.add_actor()

    def clear(self):
        for x in self.visualizers:
            x.remove_actor()

        self.visualizers = []

    def get_path_edge_tuples_set(self):
        if self.is_empty():
            return []
        return set.union(*[x.get_path_edge_tuples_set() for x in self.visualizers])


class MultiplePathPicker(MultiplePathVisualizer):
    def __init__(self, scene):
        super().__init__(scene)
        self.add_one()

    def add_one(self):
        new_one = PathPicker(self.scene, name="_picked_path_" + str(len(self.visualizers)))
        self.visualizers.append(new_one)

    def clear(self):
        super().clear()
        self.add_one()

    def set_path(self, paths, **kwargs):
        super().set_path(paths, **kwargs)

        for viz in self.visualizers:
            viz.deactivate()

    def on_start_new_path(self):
        if self.last.is_empty():
            return

        self.last.deactivate()

        self.add_one()

    def on_undo(self):
        if self.last.is_empty():
            self.visualizers.pop()
            self.last.activate()
        else:
            self.last.on_undo()

        self.scene.remove_text()
        self.scene.on_clear_result_path()

    ##############################################
    # methods where we expect just one visualizer
    ##############################################
    def get_single_point_index(self):
        if not self.is_one():
            return
        return self.last.get_single_point_index()

    def get_corresponding_edge(self):
        if not self.is_one():
            return
        return self.last.get_corresponding_edge()

    def is_intersection_point(self, idx):
        if not self.is_one():
            return
        return self.last.is_intersection_point(idx)

    ###########################################
    #  methods we simply carry on to the last
    ###########################################
    def on_close_loop(self):
        return self.last.on_close_loop()

    def on_pick_next_edge(self):
        return self.last.on_pick_next_edge()

    def on_pick_twin_edge(self):
        return self.last.on_pick_twin_edge()

    def on_pick(self, idx):
        return self.last.on_pick(idx)

