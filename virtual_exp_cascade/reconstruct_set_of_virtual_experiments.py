from samba.dcerpc.dns import soa_record


def perform(dirs, files, par, organizer):
    from settings import par, files, dirs, organizer
    from numpy import random
    import numpy as np
    from scipy.spatial import cKDTree
    from scipy.interpolate import interp1d
    import os.path
    import matplotlib.pyplot as plt

    dirpath = '/home/wgryglas/Documents/studia doktoranckie/seminaria/za2017/figures/cascade'

    from wg.tools.function import fixargument
    plt.figure = fixargument(plt.figure, override=True, figsize=(12, 6))

    modes = [m for m in organizer.load_modes()]
    modes = [modes[i] for i in range(8)]


    def argsortIndices(X, Y, geom):
        paramspace = np.array(geom.getParamList([[x,y] for x,y in zip(X, Y)], 1))
        return np.argsort(paramspace)


    base_mesh = organizer.load_mesh()
    bIds = organizer.get_mesh_boundary_ids()
    bX, bY = base_mesh.data[bIds,0], base_mesh.data[bIds,1]
    # bX, bY = organizer.load_mesh_boundary_coordinates()

    geom = organizer.load_get_geometry()
    paramSpaceBase = np.array(geom.getParamList([[x,y] for x,y in zip(bX, bY)], 1))
    sortBase = np.argsort(paramSpaceBase)

    opt_mesh_ids = organizer.load(files.optimized_mesh_ids)
    init_mesh_ids = organizer.load(files.inistial_mesh_ids)

    boundaryBaseMesh = cKDTree(np.array([bX, bY]).T)
    optXY = base_mesh.data[opt_mesh_ids]
    opt_boundary_ids = boundaryBaseMesh.query(optXY)[1]
    opt_boundary_ids_sort = argsortIndices(bX[opt_boundary_ids], bY[opt_boundary_ids], geom)

    initXY = base_mesh.data[init_mesh_ids]
    init_boundary_ids = boundaryBaseMesh.query(initXY)[1]
    init_boundary_ids = init_boundary_ids[argsortIndices(bX[init_boundary_ids], bY[init_boundary_ids], geom)]
    iParamSpace = np.array(geom.getParamList([[x,y] for x,y in initXY], 1))

    fine_mesh = organizer.load_virtual_exp_mesh()
    fine_mesh_bIds = fine_mesh.query(organizer.load_data_file(files.virtual_experiment_boundary_coords))[1]
    fine_mesh_base_ids = fine_mesh.query(np.array([bX, bY]).T)[1]
    fine_mesh_init_ids = fine_mesh.query(initXY)[1]


    fbX,fbY = fine_mesh.data[fine_mesh_bIds,0], fine_mesh.data[fine_mesh_bIds,1]
    paramSpaceFine = np.array(geom.getParamList([[x,y] for x,y in zip(fbX, fbY)], 1))
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

        x = paramSpaceBase
        y = np.mean(realizations, axis=0)
        stdderiv = np.std(realizations, axis=0)
        y1 = y - stdderiv
        y2 = y + stdderiv
        stdderiv_opt = stdderiv

        if par.save_plot_data:
            #save reconstruction on boundary
            organizer.save_plot_data(files.plot_reconstruction_result(name),
                                     {"param": x[sortBase],
                                      "x": bX[sortBase],
                                      "y_mean": y[sortBase],
                                      "y_min": y1[sortBase],
                                      "y_max": y2[sortBase],
                                      "y_exact": y_exact[sortBase],
                                      "stddev": stdderiv[sortBase]})

            organizer.save_plot_data(files.plot_reconstruction_result_in_measurement(name),
                                     {"param": x[opt_boundary_ids],
                                      "x": bX[opt_boundary_ids],
                                      "y_mean": y[opt_boundary_ids],
                                      "y_min": y1[opt_boundary_ids],
                                      "y_max": y2[opt_boundary_ids],
                                      "y_exact":y_exact[opt_boundary_ids],
                                      "stddev": stdderiv[opt_boundary_ids]})


        if par.plot:
            plt.figure()
            # for xx, yy in zip(x[opt_boundary_ids], data[fine_mesh_base_ids][opt_boundary_ids]):
            #     plt.plot([xx, xx], [yy, 1e6], "k--", lw=1)
            # plt.plot(x[opt_boundary_ids][opt_boundary_ids_sort], data[fine_mesh_base_ids][opt_boundary_ids][opt_boundary_ids_sort], 'ko', markerfacecolor='None', )

            for xx, yy in zip(x[opt_boundary_ids], data[fine_mesh_base_ids][opt_boundary_ids]):
                plt.plot([xx, xx], [-1e6, yy], "k--", lw=1, alpha=0.5)

            plt.plot(paramSpaceBase[sortBase], data[fine_mesh_base_ids][sortBase], 'k.', markersize=4, label='"Real"')

            # plt.plot(x[sortBase], y[sortBase], "k--", markerfacecolor='k', label="Reconstructed")

            plt.plot(x[sortBase], y1[sortBase], 'r--', lw=1, label='Reconstructed')
            plt.plot(x[sortBase], y2[sortBase], 'r--', lw=1)
            plt.fill_between(x[sortBase], y1[sortBase], y2[sortBase], where=y2[sortBase] >= y1[sortBase], facecolor='r', interpolate=True, alpha=0.2)

            ymin = min(y)
            ymax = max(y)
            plt.axes().set_ylim([ymin-0.1*(ymax-ymin), ymax+0.1*(ymax-ymin)])


        x = paramSpaceBase
        y = np.mean(realizations_linear, axis=0)
        stdderiv = np.std(realizations_linear, axis=0)
        y1 = y - stdderiv
        y2 = y + stdderiv

        if par.save_plot_data:
            organizer.save_plot_data(files.plot_lininterpolation_result(name),
                                     {"param": x[sortBase],
                                      "x": bX[sortBase],
                                      "y_mean": y[sortBase],
                                      "y_min": y1[sortBase],
                                      "y_max": y2[sortBase],
                                      "y_exact":y_exact[sortBase],
                                      "stddev": stdderiv[sortBase]})

            organizer.save_plot_data(files.plot_lininterpolation_result_in_measurement(name),
                                     {"param": x[init_boundary_ids],
                                      "x": bX[init_boundary_ids],
                                      "y_mean": y[init_boundary_ids],
                                      "y_min": y1[init_boundary_ids],
                                      "y_max": y2[init_boundary_ids],
                                      "y_exact": y_exact[init_boundary_ids],
                                      "stddev": stdderiv[init_boundary_ids]})

        if par.plot:
            for xx, yy in zip(x[init_boundary_ids], data[fine_mesh_base_ids][init_boundary_ids]):
                plt.plot([xx, xx], [yy, 1e6], "k--", lw=1, alpha=0.5)
            # plt.plot(x[init_boundary_ids], data[fine_mesh_base_ids][init_boundary_ids], 'ko', markerfacecolor='None', )

            # plt.plot(paramSpaceBase[sortBase], data[fine_mesh_base_ids][sortBase], 'k-', label='Viratual experiment')
            # plt.plot(x[sortBase], y[sortBase], "k--", markerfacecolor='k', label="Interpolation")

            plt.plot(x[sortBase], y1[sortBase], 'g--', lw=1, label="Interpolation")
            plt.plot(x[sortBase], y2[sortBase], 'g--', lw=1)


            #smth wrong?
            #plt.fill_between(x[sortBase], y1[sortBase], y2[sortBase], where=y2[sortBase] >= y1[sortBase], facecolor='g', interpolate=True, alpha=0.2)

            ymin = min(y)
            ymax = max(y)
            plt.axes().set_ylim([ymin-0.1*(ymax-ymin), ymax+0.1*(ymax-ymin)])

            plt.legend()
            plt.grid(True)
            plt.xlabel("Position on the profile, LE=0.5")
            plt.ylabel("Dimensionless total pressure")
            plt.show()

        if par.save_plot_data:
            organizer.save_plot_data(files.plot_exact_result(name),
                                     {"param": paramSpaceBase[sortBase],
                                      "x": bX[sortBase],
                                      "y": y_exact[sortBase]})

        if par.plot:
            plt.plot(x[sortBase], stdderiv_opt[sortBase], 'r-', lw=2, label='Reconstruction')
            plt.plot(x[sortBase], stdderiv[sortBase], 'g-', lw=2, label='Interpolation')
            plt.grid()
            plt.xlabel("Position on the profile, LE=0.5")
            plt.ylabel("Std. deviation of pressure")
            plt.legend(loc=2)
            # plt.savefig(dirpath+'/press_cascade_std_dev_'+name+'.pdf', transparent=True)
            plt.show()

