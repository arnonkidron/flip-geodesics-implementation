from scene import *


if __name__ == '__main__':
    scene = Scene('C:\\Users\\Arnon\\Desktop\\knot1.obj')

    # scene.on_pick_by_index(2067)
    # scene.on_pick_by_index(2097)
    # scene.on_pick_by_index(2094)
    # scene.path_picker.on_start_new_path()
    # scene.on_pick_by_index(2097)

    scene.on_pick_by_index(70)
    scene.on_make_geodesic()

    # scene.on_pick_by_index(2087)
    # scene.on_pick_by_index(2085)
    # scene.on_pick_by_index(2114)
    # scene.on_pick_by_index(2087)


    # scene = Scene('C:\\Users\\Arnon\\Desktop\\block.obj')
    # scene.on_pick_by_index(2087)
    # scene.on_pick_by_index(2085)
    # scene.on_pick_by_index(2114)

    # scene = Scene('C:\\Users\\Arnon\\Desktop\\horse.obj')
    # scene.on_pick_by_index(14042)
    # scene.on_pick_by_index(14257)

    # scene = Scene()
    # scene.on_pick_by_index(1733)
    # scene.on_pick_by_index(1870)
    # scene.on_pick_by_index(2078)
    # scene.on_pick_by_index(2105)

    # scene.on_flip_out()
    # scene.on_flip_out()
    # scene.on_flip_out()
    # scene.on_flip_out()
    # scene.on_flip_out()
    # scene.on_flip_out()
    # scene.on_clear()
    # scene.on_pick_by_index(2098)
    # scene.on_pick_by_index(2106)

    # scene.on_pick_by_index(1564)
    # scene.on_pick_by_index(2097)
    #
    # scene.on_make_geodesic()
    #
    # scene.on_clear()

    # scene.on_pick_by_index(1916)
    # scene.on_pick_by_index(2109)

    # DONE: check why intersection fails
    # scene.on_pick_by_index(1789)
    # scene.on_pick_by_index(1929)

    # DONE: check why intersection fails? Because of the low INTERSECTION_THRESHOLD
    # scene.on_pick_by_index(349)
    # scene.on_pick_by_index(157)
    # scene.on_flip_out()
    # scene.on_flip_out()
    # scene.on_flip_out()

    # scene.on_pick_by_index(1876)
    # scene.on_pick_by_index(1852)
    # scene.on_pick_by_index(1875)
    # scene.path_picker.on_close_loop()

    # scene.on_pick_by_index(2069)
    # scene.on_pick_by_index(2094)
    # scene.on_flip_out()
    # scene.on_flip_out()
    # scene.on_make_geodesic()


    # scene.on_pick_by_index(1741)
    # scene.on_pick_by_index(1875)
    # scene.on_make_geodesic()

    #TODO:
    # scene.on_pick_by_index(1063)
    # scene.on_pick_by_index(936)

    # scene.on_pick_by_index(14)
    # scene.on_pick_by_index(23)

    scene.show()



