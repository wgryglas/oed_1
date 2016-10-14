def perform(dirs, files, par, organizer):
    from settings import par, files, dirs, organizer
    from numpy import random
    import numpy as np
    from scipy.spatial import cKDTree
    from scipy.interpolate import interp1d

    import matplotlib.pyplot as plt

    modes = [m for m in organizer.load_modes()]

    def argsortIndices(X, Y, geom):
        paramspace = np.array(geom.getParamList([[x,y] for x,y in zip(X, Y)], 1))
        return np.argsort(paramspace)


    bIds = organizer.get_mesh_boundary_ids()
    bX, bY = organizer.load_mesh_boundary_coordinates()

    geom = organizer.load_get_geometry()
    paramSpaceBase = np.array(geom.getParamList([[x,y] for x,y in zip(bX, bY)], 1))
    sortBase = np.argsort(paramSpaceBase)

    base_mesh = organizer.load_mesh()

    opt_mesh_ids = organizer.load(files.optimized_mesh_ids)
    init_mesh_ids = organizer.load(files.inistial_mesh_ids)

    iX, iY = base_mesh.data[init_mesh_ids, 0], base_mesh.data[init_mesh_ids, 1]
    iParamSpace = np.array(geom.getParamList([[x,y] for x,y in zip(iX, iY)], 1))


    fine_mesh = organizer.load_virtual_exp_mesh()
    fine_mesh_bIds = fine_mesh.query(organizer.load_data_file(files.virtual_experiment_boundary_coords))[1]
    fine_mesh_base_ids = fine_mesh.query(np.array([bX, bY]).T)[1]
    fine_mesh_init_ids = fine_mesh.query(np.array([iX, iY]).T)[1]


    fbX,fbY = fine_mesh.data[fine_mesh_bIds,0], fine_mesh.data[fine_mesh_bIds,1]
    paramSpaceFine = np.array(geom.getParamList([[x,y] for x,y in zip(fbX, fbY)], 1))
    sortFine = np.argsort(paramSpaceFine)

    boundaryBaseMesh = cKDTree(np.array([bX, bY]).T)
    optXY = base_mesh.data[opt_mesh_ids]
    opt_boundary_ids = boundaryBaseMesh.query(optXY)[1]
    opt_boundary_ids_sort = argsortIndices(bX[opt_boundary_ids], bY[opt_boundary_ids], geom)


    initXY = base_mesh.data[init_mesh_ids]
    init_boundary_ids = boundaryBaseMesh.query(initXY)[1]
    init_boundary_ids = init_boundary_ids[argsortIndices(bX[init_boundary_ids], bY[init_boundary_ids], geom)]


    veBIds = organizer.load(files.virtual_experiment_probe_mesh_ids)

    fine_mesh_opt_ids = fine_mesh.query(optXY)[1]


    plotPntsParam = np.array(geom.getXYList(np.linspace(0, 1, 100), 1))
    boundaryMeshPlotIds = boundaryBaseMesh.query(plotPntsParam)[1]
    plotPntsParam = np.array(geom.getParamList([[x,y] for x,y in zip(bX[boundaryMeshPlotIds], bY[boundaryMeshPlotIds])], 1))
    baseMeshPlotIds = boundaryMeshPlotIds[np.argsort(plotPntsParam)]


    for datafile in files.virtual_experiment_data_files:
        reddata = organizer.load_red_file(datafile)

        data = reddata.ndata[par.optimization_variable]

        deviation = np.abs(max(data[fine_mesh_bIds]) - min(data[fine_mesh_bIds]))*par.standard_deviation/2

        bnd_data = data[fine_mesh_base_ids]

        opt_data = data[fine_mesh_opt_ids]

        realizations = list()
        realizations_linear=list()

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

            realizations_linear.append(interp1d(iParamSpace, dist_data[init_boundary_ids], kind="linear",bounds_error=False)(paramSpaceBase))

            # realizations_linear.append(dist_data[init_boundary_ids])


        plt.figure(**organizer.get_big_figure_setup())

        x = paramSpaceBase
        y = np.mean(realizations, axis=0)
        stdderiv = np.std(realizations, axis=0)
        y1 = y - stdderiv
        y2 = y + stdderiv

        # for xx, yy in zip(x[opt_boundary_ids], data[fine_mesh_base_ids][opt_boundary_ids]):
        #     plt.plot([xx, xx], [yy, 1e6], "k--", lw=1)
        plt.plot(x[opt_boundary_ids][opt_boundary_ids_sort], data[fine_mesh_base_ids][opt_boundary_ids][opt_boundary_ids_sort], 'ko', markerfacecolor='None', )

        # for xx, yy in zip(x[opt_boundary_ids], data[fine_mesh_base_ids][opt_boundary_ids]):
        #     plt.plot([xx, xx], [-1e6, yy], "k--", lw=1)

        plt.plot(paramSpaceBase[sortBase], data[fine_mesh_base_ids][sortBase], 'k-', label='Viratual experiment')
        plt.plot(x[sortBase], y[sortBase], "k--", markerfacecolor='k', label="Reconstructed")

        plt.plot(x[sortBase], y1[sortBase], 'k--', lw=1, label='Std. deviation')
        plt.plot(x[sortBase], y2[sortBase], 'k--', lw=1)
        plt.fill_between(x[sortBase], y1[sortBase], y2[sortBase], where=y2[sortBase] >= y1[sortBase], facecolor='k', interpolate=True, alpha=0.2, label='Std. deviation')

        ymin = min(y)
        ymax = max(y)
        plt.axes().set_ylim([ymin-0.1*(ymax-ymin), ymax+0.1*(ymax-ymin)])

        plt.legend()
        plt.grid(True)
        plt.savefig(dirs.figures+"/parmeter_reconstruction.png")
        plt.show()


        plt.figure(**organizer.get_big_figure_setup())

        x = paramSpaceBase
        y = np.mean(realizations_linear, axis=0)
        stdderiv = np.std(realizations_linear, axis=0)
        y1 = y - stdderiv
        y2 = y + stdderiv

        # for xx, yy in zip(x[init_boundary_ids], data[fine_mesh_base_ids][init_boundary_ids]):
        #     plt.plot([xx, xx], [yy, 1e6], "k--", lw=1)
        plt.plot(x[init_boundary_ids], data[fine_mesh_base_ids][init_boundary_ids], 'ko', markerfacecolor='None', )


        plt.plot(paramSpaceBase[sortBase], data[fine_mesh_base_ids][sortBase], 'k-', label='Viratual experiment')
        plt.plot(x[sortBase], y[sortBase], "k--", markerfacecolor='k', label="Uniform linear")

        plt.plot(x[sortBase], y1[sortBase], 'k--', lw=1, label='Std. deviation')
        plt.plot(x[sortBase], y2[sortBase], 'k--', lw=1)
        plt.fill_between(x[sortBase], y1[sortBase], y2[sortBase], where=y2[sortBase] >= y1[sortBase], facecolor='k', interpolate=True, alpha=0.2)

        ymin = min(y)
        ymax = max(y)
        plt.axes().set_ylim([ymin-0.1*(ymax-ymin), ymax+0.1*(ymax-ymin)])

        plt.legend()
        plt.grid(True)
        plt.savefig(dirs.figures+"/parmeter_linear.png")
        plt.show()
