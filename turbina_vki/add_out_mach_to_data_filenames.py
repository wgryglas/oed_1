__author__ = 'wgryglas'


def perform(dirs, files, par, organizer):
    """
    Function which addes outflow mach number to the file names
    :param dirs: object containing directories setup
    :param files: object containing files setup
    :param par: object containing parameters setup
    :return:
    """
    import glob, os
    import numpy as np
    import matplotlib.pyplot as plt
    from bearded_octo_wookie.RED.RedReader import RedTecplotFile
    from bearded_octo_wookie.RED.RedReader import GetReader
    from scipy.optimize import minimize
    from scipy.integrate import trapz
    from scipy.spatial import cKDTree

    # read geometry
    curves = GetReader(files.geom_get)
    end_curve_id = 3

    airfoil = np.array(curves.getXYList(np.linspace(0, 1, 100), 1))
    max_id = np.argmax(airfoil[:, 0])
    line_xy = np.array(curves.getXYList(np.linspace(0, 1, 100), end_curve_id))

    def dist(x):
        nLine = np.copy(line_xy)
        nLine[:, 0] -= x
        nLine -= airfoil[max_id,:]
        return min([np.dot(i, i) for i in nLine])

    dx = minimize(dist, x0=[0.], method="Nelder-Mead").x

    cord = max(airfoil[:, 0]) - min(airfoil[:, 0])

    line_xy[:, 0] = line_xy[:, 0] - dx + cord/2


    # ------- PREPARE DATA ----------------------------------------------------------------------------------------------- #

    print files.red_converged_results[0]

    init_data = RedTecplotFile(files.red_converged_results[0], useCache=True)

    mesh = cKDTree(init_data.data[:, :2])
    closest_line_ids = mesh.query(line_xy)[1]

    mesh_linexy = mesh.data[closest_line_ids, :]
    max_mesh_y = max(mesh_linexy[:, 1])
    final_ids = np.array([xy[1] < max_mesh_y for xy in line_xy], dtype=bool)
    closest_line_ids = closest_line_ids[final_ids]

    # plot extracted points
    # plt.figure()
    # plt.plot(line_xy[:, 0], line_xy[:, 1], '.')
    # plt.plot(airfoil[:, 0], airfoil[:, 1], 'r.')
    # plt.plot(airfoil[max_id, 0], airfoil[max_id, 1], 'r.', markersize=20)
    # plt.plot(mesh.data[closest_line_ids, 0], mesh.data[closest_line_ids, 1],'g.')
    # plt.show()


    class DistInterpolant:
        def __init__(self, fromXY, toXY):
            dists = np.array([(fromXY[:, 0] - xy[0])**2 + (fromXY[:, 1] - xy[1])**2 for xy in toXY])
            self.factors = np.array([d/sum(d) for d in dists])

        def interpolate(self, input_data):
            return np.array([np.dot(factor, input_data) for factor in self.factors])

    mach = DistInterpolant(mesh.data[closest_line_ids, :], line_xy)


    inputData = list()
    line_params = np.linspace(0, 1, len(line_xy))
    for i, fpath in enumerate(files.red_converged_results):
        f_data = RedTecplotFile(fpath, useCache=False)
                 # renameVariables(["X", "Y", "rho", "rhoE", "rhoU", "rhoV", "nut"]). \
                 # computeVarialbe("p", lambda x, y, rho, e, u, v, nut: (par.kappa-1.) * (e - 0.5*(u**2 + v**2)/rho)).\
                 # computeVarialbe("mach", lambda x, y, rho, e, u, v, nut, p: np.sqrt( (u/rho)**2 + (v/rho)**2) / np.sqrt(par.kappa*p/rho))

        mach_on_line = mach.interpolate(f_data.ndata.mach[closest_line_ids])
        int_avg = trapz(mach_on_line, line_params)


        # Plot output mach
        # plt.figure()
        # plt.plot(np.linspace(0, 1, len(line_xy)), mach_on_line)
        # avg_line, = plt.plot([0, 1], [int_avg, int_avg], 'g-')
        # plt.title('Out mach')
        # plt.legend([avg_line], ['averge mach %.2lf' % int_avg])
        # plt.show()

        fname = os.path.basename(fpath)
        nname = fname.split('.')[0]+("outmach%.3lf" % int_avg).replace('.', '-')+'.'+fname.split('.')[1]
        npath = os.path.dirname(fpath) + os.sep + nname
        print "renaming " + fpath+" --> " + npath
        os.rename(fpath, npath)



# if __name__ == "__main__":
#     import settings
#     perform(settings.dirs, settings.files, settings.par)