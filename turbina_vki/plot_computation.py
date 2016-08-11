
def provideData(files, par, varSelector, maxNum = None, filterFile=None):
    from bearded_octo_wookie.RED.RedReader import RedTecplotFile
    import numpy as np
    from RedBoundaryReader import readBoundaryNodes
    from scipy.spatial import cKDTree

    powKappa = (par.kappa - 1. ) /par.kappa
    p_inf = 1.

    bids = None
    bX, bY = readBoundaryNodes(files.boundary_coords)

    bX = [(min(bX) + max(bX))/2]*len(bX)
    bY = np.linspace(min(bY), max(bY), len(bX))

    count = 1
    for f in files.red_converged_results:

        if filterFile and not filterFile(f):
            continue

        if maxNum and maxNum < count:
            return

        idata = RedTecplotFile(f, useCache=True).\
                 renameVariables(["X", "Y", "rho", "rhoE", "rhoU", "rhoV", "nut"]). \
                 computeVarialbe("p", lambda x, y, rho, e, u, v, nut: (par.kappa-1.) * (e - 0.5*(u**2 + v**2)/rho)).\
                 computeVarialbe("mach", lambda x, y, rho, e, u, v, nut, p: np.sqrt( (u/rho)**2 + (v/rho)**2) / np.sqrt(par.kappa*p/rho)). \
                 computeVarialbe("mach_iso", lambda x, y, rho, e, u, v, nut, p, mach: np.sqrt(np.abs(2. / (par.kappa - 1.) * (np.abs(p_inf / p) ** (powKappa) - 1.))))

        if bids is None:
            bids = cKDTree(np.array([idata.ndata.X, idata.ndata.Y]).T).query(np.array([bX, bY]).T)[1]

        count += 1
        yield (idata.ndata.X[bids], idata.ndata.Y[bids], varSelector(idata)[bids])


def assertFile(fname, check):
    outmach = float(fname.split("outmach")[1].split(".")[0].replace("-", "."))

    ang = float(".".join(fname.split("-")[1:2]))

    if check(outmach, ang):
        print "mach="+str(outmach), "ang="+str(ang)

    return check(outmach, ang)


def provideVariableData(files, par, varSelector, maxNum, maxAng, maxMach, minMach=0, minAng=0):
    return provideData(files, par, varSelector=varSelector, maxNum=maxNum,
                       filterFile=lambda x: assertFile(x, lambda mach, ang:
                                                            minMach <= mach <= maxMach and minAng <= ang <= maxAng))


def addPlot(dirs, files, par, organizer):
    import matplotlib.pyplot as plt
    import numpy as np
    bids = organizer.get_mesh_boundary_ids()

    for idata in organizer.load_selected_red_data():
        named = idata.ndata
        x = named.X[bids]
        y = named.Y[bids]

        # x = (x*np.cos(0.46600291)-y*np.sin(0.46600291))+0.239928884974

        # x = x/ (max(x)- min(x)) - min(x)

        y = named[par.optimization_variable][bids]

        # y = y/(max(y) - min(y))

        plt.plot(x, y, ".")




def perform(dirs, files, par):
    import matplotlib.pyplot as plt

    plt.figure()

    for data in provideVariableData(files, par, varSelector=lambda d: d.ndata.mach, maxNum=20, maxMach=0.5, maxAng=95.0):
        plt.plot(data[0], data[2], ".")

    plt.show()


if __name__ == "__main__":
    from settings import *

    perform(dirs, files, par)
