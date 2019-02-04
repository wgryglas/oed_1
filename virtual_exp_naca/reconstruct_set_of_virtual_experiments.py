def perform(dirs, files, par, organizer):
    from settings import par, files, dirs, organizer
    from numpy import random
    import numpy as np
    from scipy.spatial import cKDTree
    from scipy.interpolate import interp1d
    from wg.tools.function import fixargument
    import matplotlib.pyplot as plt
    import os.path
    plt.figure = fixargument(plt.figure, override=True, figsize=(12, 6))

    modes = [m for m in organizer.load_modes()]

    b_ids = organizer.get_mesh_boundary_ids()
    bx, by = organizer.load_mesh_boundary_coordinates()

    geom = organizer.load_get_geometry()
    param_space_base = np.array(geom.getParamList([[x, y] for x, y in zip(bx, by)], 1))
    sort_base = np.argsort(param_space_base)

    base_mesh = organizer.load_mesh()

    opt_mesh_ids = organizer.load(files.optimized_mesh_ids)
    init_mesh_ids = organizer.load(files.inistial_mesh_ids)

    ix, iy = base_mesh.data[init_mesh_ids, 0], base_mesh.data[init_mesh_ids, 1]
    iParamSpace = np.array(geom.getParamList([[x, y] for x, y in zip(ix, iy)], 1))

    fine_mesh = organizer.load_virtual_exp_mesh()
    fine_mesh_b_ids = fine_mesh.query(organizer.load_data_file(files.virtual_experiment_boundary_coords))[1]
    fine_mesh_base_ids = fine_mesh.query(np.array([bx, by]).T)[1]
    fine_mesh_init_ids = fine_mesh.query(np.array([ix, iy]).T)[1]

    fbx, fby = fine_mesh.data[fine_mesh_b_ids, 0], fine_mesh.data[fine_mesh_b_ids, 1]
    param_space_fine = np.array(geom.getParamList([[x, y] for x, y in zip(fbx, fby)], 1))
    sortFine = np.argsort(param_space_fine)

    boundary_base_mesh = cKDTree(np.array([bx, by]).T)
    opt_xy = base_mesh.data[opt_mesh_ids]
    opt_boundary_ids = boundary_base_mesh.query(opt_xy)[1]

    init_xy = base_mesh.data[init_mesh_ids]
    init_boundary_ids = boundary_base_mesh.query(init_xy)[1]

    veBIds = organizer.load(files.virtual_experiment_probe_mesh_ids)

    fine_mesh_opt_ids = fine_mesh.query(opt_xy)[1]

    for datafile in files.virtual_experiment_data_files:

        name = os.path.basename(datafile).replace('.', '-')

        reddata = organizer.load_red_file(datafile)

        data = reddata.ndata[par.optimization_variable]

        deviation = np.abs(max(data[fine_mesh_b_ids]) - min(data[fine_mesh_b_ids])) * par.standard_deviation / 2

        bnd_data = data[fine_mesh_base_ids]

        opt_data = data[fine_mesh_opt_ids]

        realizations = list()
        realizations_linear = list()

        for i in range(par.virtual_exp_distortion_iteration):
            normal = random.standard_normal(len(bnd_data))
            dist_data = deviation * normal + bnd_data

            opt_dist_data = dist_data[opt_boundary_ids]

            var_modes = np.array([m.ndata[par.optimization_variable] for m in modes]).T

            A = np.array([m.ndata[par.optimization_variable][opt_mesh_ids] for m in modes]).T
            coeffs = np.linalg.lstsq(A, opt_dist_data)[0]

            result = var_modes[b_ids, :].dot(coeffs)

            realizations.append(result)

            interp = interp1d(param_space_base[init_boundary_ids], dist_data[init_boundary_ids])
            realizations_linear.append(interp(param_space_base))

        x = param_space_base
        y_exact = data[fine_mesh_base_ids][sort_base]
        y = np.mean(realizations, axis=0)
        stddev = np.std(realizations, axis=0)
        y1 = y - stddev
        y2 = y + stddev

        stddev_opt = stddev

        firstHalf = np.where(param_space_base < 0.5)
        secondHalf = np.append(np.where(param_space_base >= 0.5), 0)

        if par.save_plot_data:
            # save reconstruction on boundary
            organizer.save_plot_data(files.plot_reconstruction_result(name),
                                     {"param": x,
                                      "x": bx,
                                      "y_mean": y,
                                      "y_min": y1,
                                      "y_max": y2,
                                      "y_exact": y_exact,
                                      "stddev": stddev})

            organizer.save_plot_data(files.plot_reconstruction_result_in_measurement(name),
                                     {"param": x[opt_boundary_ids],
                                      "x": bx[opt_boundary_ids],
                                      "y_mean": y[opt_boundary_ids],
                                      "y_min": y1[opt_boundary_ids],
                                      "y_max": y2[opt_boundary_ids],
                                      "y_exact": y_exact[opt_boundary_ids],
                                      "stddev": stddev[opt_boundary_ids]})

            for half, fileName in zip([firstHalf, secondHalf], [files.plot_reconstruction_lowerCurve(name), files.plot_reconstruction_upperCurve(name)]):
                organizer.save_plot_data(fileName,
                                         {"param": x[half],
                                          "x": bx[half],
                                          "y_mean": y[half],
                                          "y_min": y1[half],
                                          "y_max": y2[half],
                                          "y_exact": y_exact[half],
                                          "stddev": stddev[half]})

        xplot = param_space_base if par.plot_in_parameter_sapce else bx

        if par.plot:
            for xx, yy in zip(xplot[opt_boundary_ids], data[fine_mesh_base_ids][opt_boundary_ids]):
                plt.plot([xx, xx], [-1e6, yy], "k--", lw=1, alpha=0.5)

            plt.plot(xplot, y1, 'r--', lw=1, label='Reconstruction')
            plt.plot(xplot, y2, 'r--', lw=1)
            for sub in [firstHalf, secondHalf]:
                plt.fill_between(xplot[sub], y1[sub], y2[sub],  where=y2[sub] >= y1[sub], facecolor='r', interpolate=True, alpha=0.2)

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
                                     {"param": x,
                                      "x": bx,
                                      "y_mean": y,
                                      "y_min": y1,
                                      "y_max": y2,
                                      "y_exact": y_exact,
                                      "stddev": stddev})

            organizer.save_plot_data(files.plot_lininterpolation_result_in_measurement(name),
                                     {"param": x[init_boundary_ids],
                                      "x": bx[init_boundary_ids],
                                      "y_mean": y[init_boundary_ids],
                                      "y_min": y1[init_boundary_ids],
                                      "y_max": y2[init_boundary_ids],
                                      "y_exact": y_exact[init_boundary_ids],
                                      "stddev": stddev[init_boundary_ids]})

            for half, fileName in zip([firstHalf, secondHalf], [files.plot_lininterpolation_lowerCurve(name), files.plot_lininterpolation_upperCurve(name)]):
                organizer.save_plot_data(fileName,
                                         {"param": x[half],
                                          "x": bx[half],
                                          "y_mean": y[half],
                                          "y_min": y1[half],
                                          "y_max": y2[half],
                                          "y_exact": y_exact[half],
                                          "stddev": stddev[half]})

        if par.plot:
            for xx, yy in zip(xplot[init_boundary_ids], data[fine_mesh_base_ids][init_boundary_ids]):
                plt.plot([xx, xx], [yy, 1e6], "k--", lw=1, alpha=0.5)

            plt.plot(xplot, y1, 'g--', lw=1, label='Interpolation')
            plt.plot(xplot, y2, 'g--', lw=1)

            for sub in [firstHalf, secondHalf]:
                plt.fill_between(xplot[sub], y1[sub], y2[sub],  where=y2[sub] >= y1[sub], facecolor='g', interpolate=True, alpha=0.2)
            #plt.fill_between(xplot, y1, y2, where=y2 >= y1, facecolor='g', interpolate=True, alpha=0.2)

            ymin = min(y)
            ymax = max(y)
            plt.axes().set_ylim([ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin)])

        if par.save_plot_data:
            organizer.save_plot_data(files.plot_exact_result(name),
                                     {"param": param_space_base[sort_base],
                                      "x": bx,
                                      "y": data[fine_mesh_base_ids][sort_base]})

        if par.plot:
            plt.plot(xplot[sort_base], data[fine_mesh_base_ids][sort_base], '.', lw=1, color="black",
                     label='"Real"')

            plt.legend(loc=2 if par.plot_in_parameter_sapce else 4)

            plt.xlabel("Position on the profile, LE=0.5" if par.plot_in_parameter_sapce else "X coordinate")
            plt.ylabel(par.optimization_variable)
            # plt.savefig(dirpath+'/press_profile_airfoil_'+name+'.pdf', transparent=True)
            plt.xlim([-0.1, 1.1])
            plt.show()

        if par.plot:
            plt.plot(xplot, stddev_opt / par.standard_deviation, 'ro-', lw=2, label='Reconstruction')
            plt.plot(xplot, stddev / par.standard_deviation, 'gv-', lw=2, label='Interpolation')
            plt.grid()
            plt.xlabel("Position on the profile, LE=0.5" if par.plot_in_parameter_sapce else "X coordinate")
            plt.ylabel("Std. deviation of pressure")
            plt.legend(loc=2 if par.plot_in_parameter_sapce else 3)
            # plt.savefig(dirpath+'/press_airfoil_std_dev_'+name+'.pdf', transparent=True)
            plt.xlim([-0.1, 1.1])
            plt.show()
