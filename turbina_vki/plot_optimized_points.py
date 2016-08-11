__author__ = 'wgryglas'


def perform(dirs, files, par, organizer):

    import matplotlib.pyplot as plt

    mesh = organizer.load_mesh()
    opt_mesh_ids = organizer.load(files.optimized_mesh_ids)
    probe_mesh_ids = organizer.load(files.probe_mesh_ids)
    start_mesh_ids = organizer.load(files.inistial_mesh_ids)
    bx, by = organizer.load_mesh_boundary_coordinates()
    sort = organizer.sort_points(bx, by)
    bx = bx[sort]
    by = by[sort]

    plt.figure()
    plt.plot(bx, by, '-', color="gray", markersize=3)
    plt.plot(mesh.data[opt_mesh_ids, 0], mesh.data[opt_mesh_ids, 1], 'r.', markersize=35)
    plt.plot(mesh.data[start_mesh_ids, 0], mesh.data[start_mesh_ids, 1], '.', color="blue", markersize=25)
    plt.plot(mesh.data[probe_mesh_ids, 0], mesh.data[probe_mesh_ids, 1], '.', color="black", markersize=15)
    plt.show()


if __name__ == "__main__":
    from settings import *
    perform(dirs, files, par, organizer)