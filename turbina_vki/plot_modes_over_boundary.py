

def perform(dirs, files, par, organizer):
    import matplotlib.pyplot as plt

    x, y = organizer.load_mesh_boundary_coordinates()
    sort = organizer.sort_points(x,y)
    x = x[sort]
    y = y[sort]

    bids = organizer.get_mesh_boundary_ids()[sort]

    for i, m in enumerate(organizer.load_modes()):
        plt.plot(x, m.ndata.mach_iso[bids], "-", color="black", markersize=3)
        plt.title(r'$mode\,'+str(i)+'$')
        plt.xlabel(r'$\frac{x}{chord\, length} [-]$')
        plt.ylabel(r'$Mode\, value [-]$')
        plt.axes().grid(True)
        plt.show()



if __name__ == "__main__":
    from settings import *
    perform(dirs, files, par, organizer)




