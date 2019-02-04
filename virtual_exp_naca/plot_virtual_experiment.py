def perform(dirs, files, par, organizer):
    from bearded_octo_wookie.RED import POD
    import matplotlib.pyplot as plt
    import numpy as np

    bcoords = np.array( organizer.load_experiment_points() ).T

    ve_mesh = organizer.load_virtual_exp_mesh()
    bids = ve_mesh.query(bcoords)[1]
    sort = organizer.sort_points(bcoords[:, 0], bcoords[:, 1], organizer.load_get_geometry())
    bids = bids[sort]

    #plt.plot(ve_mesh.data[bids, 0], ve_mesh.data[bids, 1], '.')
    provider = organizer.load_virtual_exp_data()
    data = provider.next()

    print data.variables

    f, (a1, a2) = plt.subplots(1, 2)

    ref_point = np.argmin(data.ndata.X)
    print "Inlet pressure", data.ndata.p[ref_point]
    print "Inlet cp", data.ndata.cp[ref_point]
    print "Inlet dyn. pressure", (data.ndata.rhoU[ref_point]**2 + data.ndata.rhoV[ref_point]**2) / data.ndata.rho[ref_point] / 2

    a1.plot(ve_mesh.data[bids, 0], data.ndata.cp[bids])
    a2.plot(ve_mesh.data[bids, 0], data.ndata.p[bids])

    #plt.ylim([1, -2])
    plt.show()

if __name__ == "__main__":
    from settings import *
    perform(dirs, files, par, organizer)