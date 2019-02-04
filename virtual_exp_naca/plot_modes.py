def perform(dirs, files, par, organizer):
    from bearded_octo_wookie.RED import POD
    import matplotlib.pyplot as plt

    bids = organizer.get_mesh_boundary_ids()

    #mode = organizer.load_mode(1)
    for mode in organizer.load_modes():
        d = mode.ndata
        plt.plot(d.X[bids], d.cp[bids])
    plt.show()

if __name__ == "__main__":
    from settings import *
    perform(dirs, files, par, organizer)
