
def perform(dirs, files, par, organizer):
    import matplotlib.pyplot as plt
    from wg.tools.plot import latex_xticks

    boundary_ids = organizer.get_mesh_boundary_ids()
    bX, bY = organizer.load_mesh_boundary_coordinates()
    bsort = organizer.sort_points(bX, bY)


    boundary_data = {}
    for opt_var, num_modes, pressure, tecplot in organizer.load_reconstruction():
        boundary_data[pressure] = tecplot.ndata[par.optimization_variable][boundary_ids]

    pX, pY = organizer.tranform_to_mesh_csys(*organizer.load_experiment_points())
    exp_data = {}
    for p, data in organizer.load_selected_experiment_data():
        exp_data[p] = data[par.optimization_variable]

    opt_pnts = organizer.load(files.optimized_probe_ids)


    mapping = organizer.load_pressure_to_mach_mapping()

    plt.rc('text', usetex=True)

    for p in exp_data:
        plt.figure()
        l1, = plt.plot(bX[bsort], boundary_data[p][bsort], "-r", markersize=3)
        l2, = plt.plot(pX[opt_pnts], exp_data[p][opt_pnts], ".", color="black", markersize=30)
        l3, = plt.plot(pX, exp_data[p], ".", color="green", markersize=15)

        plt.title(r'$Outflow\, mach\, number\,'+str(mapping[p])+'$')
        plt.xlabel(r'$\frac{x}{chord\, length} [-]$')
        plt.ylabel(r'$Mach [-]$')
        plt.axes().grid(True)
        #plt.legend([l1, l2, l3], [r'$Reconstructed\, data$', r'Experiment\, data\, used\, in\, reconstruction$', r'$All\, measured\, data$'])

        plt.show()



    # for i, (rec, exp_data, value) in enumerate(zip(reconstructs, exp_datas, tankPresures)):
    #     plt.figure()
    #     #col = cmap(float(i)/len(exp_datas))
    #     col = [1., 0., 0., 1.]
    #
    #     #plt.plot(pX2[opt_pnts], exp_data[probe_sort][opt_pnts], '.', color="yellow", markersize=20)
    #     #plt.plot(pX2, exp_data[probe_sort], '.', color=shift_color(col, 0.25), markersize=15)
    #
    #     plt.plot(bX, rec[boundary_ids][bsort], '-', color="red", markersize=2)
    #
    #     plt.plot(pX2[opt_pnts], exp_data[opt_pnts], '.', color="black", markersize=20)
    #
    #     plt.plot(pX[probe_sort], exp_data[probe_sort], '.', color="green", markersize=10)
    #
    #     plt.title("Pressure in tank "+str(value))
    #     # plt.ylim([0, 1.2])
    #     plt.show()

if __name__ == "__main__":
    from settings import *
    perform(dirs, files, par, organizer)

