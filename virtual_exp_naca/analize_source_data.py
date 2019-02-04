

def perform(dirs, files, par, organizer):
    from bearded_octo_wookie.RED import POD
    from wg.tools.system import clear_dir
    import matplotlib.pyplot as plt

    # par.useCache = False

    alldata = organizer.load_selected_red_data()

    for i, d in enumerate(alldata):
        data = d.ndata

        print data.variables

        dynamicPressure = (data.rhoU**2 + data.rhoV**2) / data.rho / 2
        #dynamicPressure = data.rhoU**2 + data.rhoV**2 / 2

        # print data.rhoU[0] / data.rho[0]
        # print data.rhoU[0] / data.rho[0] / np.sqrt(par.kappa * p[0] / data.rho[0])

        ids = organizer.get_mesh_boundary_ids()

        # plt.plot(data.X, data.Y, 'k.', markersize=1)
        # plt.plot(data.X[ids], data.Y[ids], 'r-')
        # plt.plot(data.X[0], data.Y[0], 'g.', markersize=3)

        inVel = np.sqrt(data.rhoU[0]**2 + data.rhoV[0]**2) / data.rho[0]

        if inVel < 0.699:
            print "Inlet velocity", inVel
            print "Inlet density", data.rho[0]
            print "Inlet pressure", data.p[0]
            print "-------------"
            plt.plot(data.X[ids], data.cp[ids])

        if i > 10: break


    plt.ylim([1, -2])
    plt.show()




if __name__ == "__main__":
    from settings import *
    perform(dirs, files, par, organizer)