def perform(dirs, files, par, organizer):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.spatial import cKDTree
    from wg.tools.array import spread_id

    # ------- PREPARE DATA ----------------------------------------------------------------------------------------------- #
    np.seterr(all='raise')

    # read mesh and necessery data
    mesh = organizer.load_mesh()

    boundary_ids = organizer.get_mesh_boundary_ids()
    bX = mesh.data[boundary_ids, 0]
    bY = mesh.data[boundary_ids, 1]


    # read probe points
    pX, pY = organizer.get_probe_positions()


    # find closest indices of boundary nodes to the probe coords
    bProbeIds = cKDTree(np.array([bX, bY]).T).query(np.array([pX, pY]).T)[1]
    probe_mesh_ids = boundary_ids[bProbeIds]
    pX2 = mesh.data[probe_mesh_ids, 0]
    pY2 = mesh.data[probe_mesh_ids, 1]


    # find boundary and probe nodes parametrization (for setting appropriate order)
    bsort = organizer.sort_points(bX, bY)

    # find sorting for probes according to curve parametrization
    probe_sort = organizer.sort_points(pX2, pY2)

    def getSubFisherMatrix(modes, pntIds):
        X = np.matrix(modes[pntIds, :])
        return X.T*X

    def computeACriterion(M):
        try:
            D = np.linalg.eigvalsh(M)
            tr = np.ravel(np.sum(1./D))
            if tr < 0:
                tr = np.inf
        except Exception as e:
            print e
            tr = np.inf
        return tr

    # Perform optimization

    # Take initial positions:
    # find closest probe point to equidistributed points
    from scipy.spatial import cKDTree
    probeMesh = cKDTree(np.array([pX2, pY2]).T)
    equalProbeIds = probeMesh.query(organizer.get_equaly_distributed_points_over_boundary(par.num_measure_pnts))[1]

    start_probe_ids = equalProbeIds
    start_mesh_ids = probe_mesh_ids[start_probe_ids]

    organizer.save(files.inistial_mesh_ids, start_mesh_ids)
    organizer.save(files.inistial_probe_ids, start_probe_ids)


    opt_pnts = np.copy(start_probe_ids)
    n_pnts = len(start_probe_ids)


    # Load modes
    modes = [organizer.load_mode(i) for i in range(par.num_modes)]

    maxTest = 5
    mintr = float('inf')

    # Calc appropriate input data
    varmodes = np.array([m.ndata[par.optimization_variable] for m in modes]).T

    for loopId in range(maxTest):

        moved = False

        for i in range(n_pnts):
            mask = np.ones(n_pnts, dtype=bool)
            mask[i] = False
            tmp_opt = opt_pnts[mask]

            locmin = (-1, mintr)

            print 'checking', i, 'point with boundary position ', opt_pnts[i], 'and tr=', mintr

            traces = list()

            for j, _ in enumerate(probe_mesh_ids):
                if j in tmp_opt:
                    continue

                tmp_opt2 = np.append(tmp_opt, j)

                testPos = probe_mesh_ids[tmp_opt2]
                criterion = computeACriterion(getSubFisherMatrix(varmodes, testPos))

                if criterion < 0:
                    continue

                if criterion < locmin[1]:
                    locmin = (j, criterion)

            if locmin[0] != -1:
                opt_pnts[i] = locmin[0]
                mintr = locmin[1]
                moved = True
                print 'moved to', locmin[0], 'tr=', locmin[1]

        print '---- END OF LOOP', loopId, '-----------------------------------------------------'

        if not moved:
            break


    opt_mesh_ids = probe_mesh_ids[opt_pnts]

    organizer.save(files.optimized_probe_ids, opt_pnts)
    organizer.save(files.optimized_mesh_ids, opt_mesh_ids)
    organizer.save(files.probe_mesh_ids, probe_mesh_ids)


    # Plot result of optimization
    # plt.figure()
    # plt.plot(mesh.data[opt_mesh_ids, 0], mesh.data[opt_mesh_ids, 1], 'r.', markersize=30)
    # plt.plot(mesh.data[start_mesh_ids, 0], mesh.data[start_mesh_ids, 1], 'g.', markersize=20)
    # plt.plot(mesh.data[probe_mesh_ids, 0], mesh.data[probe_mesh_ids, 1], 'b.')
    # plt.show()


if __name__ == "__main__":
    from settings import *
    perform(dirs, files, par, organizer)