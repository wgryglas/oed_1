
def perform(dirs, files, par, organizer):
    import numpy as np
    from scipy.spatial import cKDTree
    reddata = organizer.load_red_file(files.virtual_experiment_data_files[0])

    data = reddata.ndata

    geom = organizer.load_get_geometry()

    listXY = np.array(geom.getXYList(np.linspace(0, 1, 10e3), 1))

    mesh = cKDTree(np.array([data.X, data.Y]).T)

    nodes = mesh.query(listXY)[1]
    nodes = np.unique(nodes)

    print "Identified", len(nodes), "nodes"

    import matplotlib.pyplot as plt
    plt.plot(data.X[nodes], data.Y[nodes], ".")
    plt.show()

    np.savetxt(files.virtual_experiment_boundary_coords, np.array([data.X[nodes], data.Y[nodes]]).T)





