def perform(dirs, files, par, organizer):
    from settings import par, files, dirs, organizer
    from numpy import random
    import numpy as np
    from scipy.spatial import cKDTree

    import matplotlib.pyplot as plt

    modes = [m for m in organizer.load_modes()]

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

    initXY = base_mesh.data[init_mesh_ids]
    init_boundary_ids = boundaryBaseMesh.query(initXY)[1]


    veBIds = organizer.load(files.virtual_experiment_probe_mesh_ids)

    fine_mesh_opt_ids = fine_mesh.query(optXY)[1]

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

            realizations_linear.append(dist_data[init_boundary_ids])

            # for c, m in zip(coeffs, modes):
            #     result += c*m.ndata[par.optimization_variable][bIds]

            # plt.plot(paramSpaceBase[sortBase], result[sortBase], '-', color=(0.5,0.5,0.5,0.1))

            # plt.plot(reddata.ndata.X[fine_mesh_opt_ids], opt_data, '.', color="red", markersize=10)
            # plt.plot(reddata.ndata.X[fine_mesh_opt_ids], dist_data, '.', color="pink", markersize=5)

            # plt.plot(bX, data[veBIds], '.', color="black", markersize=5)

        x = paramSpaceBase
        y = np.mean(realizations, axis=0)
        stdderiv = np.std(realizations, axis=0)
        y1 = y - stdderiv
        y2 = y + stdderiv

        # for xx, yy in zip(x[init_boundary_ids], data[fine_mesh_base_ids][init_boundary_ids]):
        #     plt.plot([xx, xx], [yy, 1e6], "k--", lw=1)

        for xx, yy in zip(x[opt_boundary_ids], data[fine_mesh_base_ids][opt_boundary_ids]):
            plt.plot([xx, xx], [-1e6, yy], "k--", lw=1)

        plt.plot(paramSpaceBase[sortBase], data[fine_mesh_base_ids][sortBase],'-', color="black", label='"Real"')
        plt.plot(x, y, ".", color="blue")

        plt.plot(x[sortBase], y1[sortBase], 'b--', lw=1)
        plt.plot(x[sortBase], y2[sortBase], 'b--', lw=1)
        plt.fill_between(x[sortBase], y1[sortBase], y2[sortBase], where=y2[sortBase] >= y1[sortBase], facecolor='b', interpolate=True, alpha=0.2)

        ymin = min(y)
        ymax = max(y)
        plt.axes().set_ylim([ymin-0.1*(ymax-ymin), ymax+0.1*(ymax-ymin)])

        plt.show()

        x = paramSpaceBase[init_boundary_ids]
        y = np.mean(realizations_linear, axis=0)
        stdderiv = np.std(realizations_linear, axis=0)
        y1 = y - stdderiv
        y2 = y + stdderiv

        for xx, yy in zip(x, data[fine_mesh_base_ids][init_boundary_ids]):
            plt.plot([xx, xx], [yy, 1e6], "k--", lw=1)

        # for xx, yy in zip(x[opt_boundary_ids], data[fine_mesh_base_ids][opt_boundary_ids]):
        #     plt.plot([xx, xx], [-1e6, yy], "k--", lw=1)

        plt.plot(paramSpaceBase[sortBase], data[fine_mesh_base_ids][sortBase],'-', color="black", label='"Real"')
        plt.plot(x, y, ".", color="blue")

        plt.plot(x, y1, 'b--', lw=1)
        plt.plot(x, y2, 'b--', lw=1)
        plt.fill_between(x, y1, y2, where=y2 >= y1, facecolor='b', interpolate=True, alpha=0.2)

        ymin = min(y)
        ymax = max(y)
        plt.axes().set_ylim([ymin-0.1*(ymax-ymin), ymax+0.1*(ymax-ymin)])

        plt.show()

        # x = bX
        # plt.plot(reddata.ndata.X[fine_mesh_base_ids][sortBase], data[fine_mesh_base_ids][sortBase],'b-', label='"Real"')
        # plt.plot(x, y, ".", color="black")
        # plt.plot(x, y1 , 'b--', lw=1)
        # plt.plot(x, y2 , 'b--', lw=1)
        # plt.fill_between(x, y1, y2, where=y2>=y1, facecolor='b', interpolate=True, alpha=0.2)
        # plt.show()