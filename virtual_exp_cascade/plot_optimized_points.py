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

    plt.figure(figsize=(10, 5))
    plt.plot(bx, by, '-', color="blue", markersize=3)
    # plt.plot(mesh.data[probe_mesh_ids, 0], mesh.data[probe_mesh_ids, 1], '.', color="black", markersize=10)
    # plt.plot(mesh.data[opt_mesh_ids, 0], mesh.data[opt_mesh_ids, 1], '.', color="black", markersize=25)
    plt.plot(mesh.data[opt_mesh_ids, 0], mesh.data[opt_mesh_ids, 1], 'o', color="green", markersize=15)
    plt.plot(mesh.data[start_mesh_ids, 0], mesh.data[start_mesh_ids, 1], 'rx', markersize=15)
    plt.axes().set_aspect('equal', 'datalim')
    plt.axes().set_xlim([-0.05, 1.05])
    plt.grid(True)
    dirpath = '/home/wgryglas/Documents/konferencje/rok2016/porto/presentation/img'
    plt.savefig(dirpath+'/cascade_opt_pos.pdf', transparent=True)
    plt.show()



if __name__ == "__main__":
    from settings import *
    perform(dirs, files, par, organizer)