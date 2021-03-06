from samba.dcerpc.dns import soa_record


def perform(dirs, files, par, organizer):
    from settings import par, files, dirs, organizer
    from numpy import random
    import numpy as np
    from scipy.spatial import cKDTree
    from scipy.interpolate import interp1d
    import os.path
    import matplotlib.pyplot as plt

    from wg.tools.function import fixargument
    plt.figure = fixargument(plt.figure, override=True, figsize=(12, 6))

    modes = [m for m in organizer.load_modes()]
    modes = [modes[i] for i in range(8)]


    def argsortIndices(X, Y, geom):
        paramspace = np.array(geom.getParamList([[x, y] for x, y in zip(X, Y)], 1))
        return np.argsort(paramspace)


    base_mesh = organizer.load_mesh()
    bIds = organizer.get_mesh_boundary_ids()
    bX, bY = base_mesh.data[bIds,0], base_mesh.data[bIds,1]
    # bX, bY = organizer.load_mesh_boundary_coordinates()

    geom = organizer.load_get_geometry()
    paramSpaceBase = np.array(geom.getParamList([[x, y] for x, y in zip(bX, bY)], 1))
    sortBase = np.argsort(paramSpaceBase)

    opt_mesh_ids = organizer.load(files.optimized_mesh_ids)
    init_mesh_ids = organizer.load(files.inistial_mesh_ids)

    boundaryBaseMesh = cKDTree(np.array([bX, bY]).T)
    optXY = base_mesh.data[opt_mesh_ids]
    opt_boundary_ids = boundaryBaseMesh.query(optXY)[1]
    opt_boundary_ids_sort = opt_boundary_ids[argsortIndices(bX[opt_boundary_ids], bY[opt_boundary_ids], geom)]

    initXY = base_mesh.data[init_mesh_ids]
    init_boundary_ids = boundaryBaseMesh.query(initXY)[1]
    init_boundary_ids_sort = init_boundary_ids[argsortIndices(bX[init_boundary_ids], bY[init_boundary_ids], geom)]
    iParamSpace = np.array(geom.getParamList([[x, y] for x, y in initXY], 1))

    fine_mesh = organizer.load_virtual_exp_mesh()
    fine_mesh_bIds = fine_mesh.query(organizer.load_data_file(files.virtual_experiment_boundary_coords))[1]
    fine_mesh_base_ids = fine_mesh.query(np.array([bX, bY]).T)[1]
    fine_mesh_init_ids = fine_mesh.query(initXY)[1]


    fbX,fbY = fine_mesh.data[fine_mesh_bIds,0], fine_mesh.data[fine_mesh_bIds,1]
    paramSpaceFine = np.array(geom.getParamList([[x, y] for x, y in zip(fbX, fbY)], 1))
    sortFine = np.argsort(paramSpaceFine)



    veBIds = organizer.load(files.virtual_experiment_probe_mesh_ids)

    fine_mesh_opt_ids = fine_mesh.query(optXY)[1]


    plotPntsParam = np.array(geom.getXYList(np.linspace(0, 1, 100), 1))
    boundaryMeshPlotIds = boundaryBaseMesh.query(plotPntsParam)[1]
    plotPntsParam = np.array(geom.getParamList([[x, y] for x, y in zip(bX[boundaryMeshPlotIds], bY[boundaryMeshPlotIds])], 1))
    baseMeshPlotIds = boundaryMeshPlotIds[np.argsort(plotPntsParam)]


    for datafile in files.virtual_experiment_data_files:

        name = os.path.basename(datafile).replace('.', '-')

        reddata = organizer.load_red_file(datafile)

        data = reddata.ndata[par.optimization_variable]

        deviation = np.abs(max(data[fine_mesh_bIds]) - min(data[fine_mesh_bIds]))*par.standard_deviation/2

        bnd_data = data[fine_mesh_base_ids]

        opt_data = data[fine_mesh_opt_ids]

        realizations = list()
        realizations_linear = list()

        for i in range(par.virtual_exp_distortion_iteration):
            normal = random.standard_normal(len(bnd_data))
            # normal = random.standard_normal(len(opt_data))
            dist_data = deviation * normal + bnd_data

            opt_dist_data = dist_data[opt_boundary_ids]

            varmodes = np.array([m.ndata[par.optimization_variable] for m in modes]).T

            A = np.array([m.ndata[par.optimization_variable][opt_mesh_ids] for m in modes]).T
            coeffs = np.linalg.lstsq(A, opt_dist_data)[0]

            result = varmodes[bIds, :].dot(coeffs)

            realizations.append(result)

            interp = interp1d(paramSpaceBase[init_boundary_ids], dist_data[init_boundary_ids], bounds_error=False)
            realizations_linear.append(interp(paramSpaceBase))
            # realizations_linear.append(interp1d(iParamSpace, dist_data[init_boundary_ids], kind="linear", bounds_error=False)(paramSpaceBase))

            # realizations_linear.append(dist_data[init_boundary_ids])

        y_exact = data[fine_mesh_base_ids]
        y = np.mean(realizations, axis=0)
        stddev = np.std(realizations, axis=0)
        y1 = y - stddev
        y2 = y + stddev

        stddev_opt = stddev

        firstHalf = np.where(paramSpaceBase[sortBase] < 0.5)
        secondHalf = np.append(np.where(paramSpaceBase[sortBase] >= 0.5), 0)

        if par.save_plot_data:
            # save reconstruction on boundary
            organizer.save_plot_data(files.plot_reconstruction_result(name),
                                     {"param": paramSpaceBase[sortBase],
                                      "x": bX[sortBase],
                                      "y_mean": y[sortBase],
                                      "y_min": y1[sortBase],
                                      "y_max": y2[sortBase],
                                      "y_exact": y_exact[sortBase],
                                      "stddev": stddev[sortBase]})

            organizer.save_plot_data(files.plot_reconstruction_result_in_measurement(name),
                                     {"param": paramSpaceBase[opt_boundary_ids_sort],
                                      "x": bX[opt_boundary_ids_sort],
                                      "y_mean": y[opt_boundary_ids_sort],
                                      "y_min": y1[opt_boundary_ids_sort],
                                      "y_max": y2[opt_boundary_ids_sort],
                                      "y_exact": y_exact[opt_boundary_ids_sort],
                                      "stddev": stddev[opt_boundary_ids_sort]})

            for half, fileName in zip([firstHalf, secondHalf], [files.plot_reconstruction_lowerCurve(name), files.plot_reconstruction_upperCurve(name)]):
                organizer.save_plot_data(fileName,
                                         {"param": paramSpaceBase[sortBase][half],
                                          "x": bX[sortBase][half],
                                          "y_mean": y[sortBase][half],
                                          "y_min": y1[sortBase][half],
                                          "y_max": y2[sortBase][half],
                                          "y_exact": y_exact[sortBase][half],
                                          "stddev": stddev[sortBase][half]})

        xplot = paramSpaceBase if par.plot_in_parameter_space else bX

        if par.plot:
            for xx, yy in zip(xplot[opt_boundary_ids_sort], data[fine_mesh_base_ids][opt_boundary_ids_sort]):
                plt.plot([xx, xx], [-1e6, yy], "k--", lw=1, alpha=0.5)

            plt.plot(xplot[sortBase], y1[sortBase], 'r--', lw=1, label='Reconstruction')
            plt.plot(xplot[sortBase], y2[sortBase], 'r--', lw=1)
            # for sub in [firstHalf, secondHalf]:
            #     plt.fill_between(xplot[sub], y1[sub], y2[sub],  where=y2[sub] >= y1[sub], facecolor='r', interpolate=True, alpha=0.2)

            ymin = min(y)
            ymax = max(y)
            plt.axes().set_ylim([ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin)])

        y = np.mean(realizations_linear, axis=0)
        stddev = np.std(realizations_linear, axis=0)
        y1 = y - stddev
        y2 = y + stddev

        if par.save_plot_data:
            # save interpolation on boundary
            organizer.save_plot_data(files.plot_lininterpolation_result(name),
                                     {"param": paramSpaceBase[sortBase],
                                      "x": bX[sortBase],
                                      "y_mean": y[sortBase],
                                      "y_min": y1[sortBase],
                                      "y_max": y2[sortBase],
                                      "y_exact": y_exact[sortBase],
                                      "stddev": stddev[sortBase]})

            organizer.save_plot_data(files.plot_lininterpolation_result_in_measurement(name),
                                     {"param": paramSpaceBase[init_boundary_ids_sort],
                                      "x": bX[init_boundary_ids_sort],
                                      "y_mean": y[init_boundary_ids_sort],
                                      "y_min": y1[init_boundary_ids_sort],
                                      "y_max": y2[init_boundary_ids_sort],
                                      "y_exact": y_exact[init_boundary_ids_sort],
                                      "stddev": stddev[init_boundary_ids_sort]})

            for half, fileName in zip([firstHalf, secondHalf], [files.plot_lininterpolation_lowerCurve(name), files.plot_lininterpolation_upperCurve(name)]):
                organizer.save_plot_data(fileName,
                                         {"param": paramSpaceBase[half],
                                          "x": bX[half],
                                          "y_mean": y[half],
                                          "y_min": y1[half],
                                          "y_max": y2[half],
                                          "y_exact": y_exact[half],
                                          "stddev": stddev[half]})

        if par.plot:
            for xx, yy in zip(xplot[init_boundary_ids], data[fine_mesh_base_ids][init_boundary_ids]):
                plt.plot([xx, xx], [yy, 1e6], "k--", lw=1, alpha=0.5)

            plt.plot(xplot[sortBase], y1[sortBase], 'g--', lw=1, label='Interpolation')
            plt.plot(xplot[sortBase], y2[sortBase], 'g--', lw=1)

            # for sub in [firstHalf, secondHalf]:
            #     plt.fill_between(xplot[sub], y1[sub], y2[sub],  where=y2[sub] >= y1[sub], facecolor='g', interpolate=True, alpha=0.2)

            #plt.fill_between(xplot, y1, y2, where=y2 >= y1, facecolor='g', interpolate=True, alpha=0.2)

            ymin = min(y)
            ymax = max(y)
            plt.axes().set_ylim([ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin)])

        if par.save_plot_data:
            organizer.save_plot_data(files.plot_exact_result(name),
                                     {"param": paramSpaceBase[sortBase],
                                      "x": bX[sortBase],
                                      "y": data[fine_mesh_base_ids][sortBase]})

        if par.plot:
            plt.plot(xplot[sortBase], data[fine_mesh_base_ids][sortBase], '.', lw=1, color="black", label='"Real"')

            plt.legend(loc=2 if par.plot_in_parameter_space else 4)

            plt.xlabel("Position on the profile, LE=0.5" if par.plot_in_parameter_space else "X coordinate")
            plt.ylabel(par.optimization_variable)
            # plt.savefig(dirpath+'/press_profile_airfoil_'+name+'.pdf', transparent=True)
            plt.xlim([-0.1, 1.1])
            plt.show()

        if par.plot:
            plt.plot(xplot[sortBase], (stddev_opt / par.standard_deviation)[sortBase], 'ro-', lw=2, label='Reconstruction')
            plt.plot(xplot[sortBase], (stddev / par.standard_deviation)[sortBase], 'gv-', lw=2, label='Interpolation')
            plt.grid()
            plt.xlabel("Position on the profile, LE=0.5" if par.plot_in_parameter_space else "X coordinate")
            plt.ylabel("Std. deviation of pressure")
            plt.legend(loc=2 if par.plot_in_parameter_space else 3)
            # plt.savefig(dirpath+'/press_airfoil_std_dev_'+name+'.pdf', transparent=True)
            plt.xlim([-0.1, 1.1])
            plt.show()