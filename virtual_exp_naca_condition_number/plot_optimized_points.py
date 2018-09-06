__author__ = 'wgryglas'


def perform(dirs, files, par, organizer):

    import matplotlib.pyplot as plt
    import pandas as pd

    mesh = organizer.load_mesh()
    opt_mesh_ids = organizer.load(files.optimized_mesh_ids)
    probe_mesh_ids = organizer.load(files.probe_mesh_ids)
    start_mesh_ids = organizer.load(files.inistial_mesh_ids)
    bx, by = organizer.load_mesh_boundary_coordinates()
    sort = organizer.sort_points(bx, by)
    bx = bx[sort]
    by = by[sort]

    if par.plot:
        plt.figure(figsize=(10, 2.5))
        plt.plot(bx, by, '-', color="blue", markersize=3)
        # plt.plot(mesh.data[probe_mesh_ids, 0], mesh.data[probe_mesh_ids, 1], '.', color="black", markersize=10)
        # plt.plot(mesh.data[opt_mesh_ids, 0], mesh.data[opt_mesh_ids, 1], '.', color="black", markersize=25)
        plt.plot(mesh.data[opt_mesh_ids, 0], mesh.data[opt_mesh_ids, 1], 'o', color="green", markersize=15)
        plt.plot(mesh.data[start_mesh_ids, 0], mesh.data[start_mesh_ids, 1], 'rx', markersize=15)
        plt.axes().set_xlim([-0.05, 1.05])
        plt.axes().set_aspect('equal', 'datalim')
        # plt.axes().set_ylim([-0.06, 0.06])
        plt.grid(True)

        #plt.savefig(dirpath+'/naca_opt_pnts.pdf', transparent=True)
        plt.show()


    if par.save_plot_data:
        organizer.save_plot_data(files.plot_profile, {"x": bx, "y": by})
        organizer.save_plot_data(files.plot_init_opt_pnts, {"xopt": mesh.data[opt_mesh_ids,   0],
                                                            "yopt": mesh.data[opt_mesh_ids,   1],
                                                            "xeq":  mesh.data[start_mesh_ids, 0],
                                                            "yeq":  mesh.data[start_mesh_ids, 1]})



if __name__ == "__main__":
    from settings import *
    perform(dirs, files, par, organizer)