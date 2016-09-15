

def perform(dirs, files, par, organizer):
    import numpy as np
    modes = [m for m in organizer.load_modes()]

    opt_mesh_ids = organizer.load(files.optimized_mesh_ids)
    opt_probe_pnts = organizer.load(files.optimized_probe_ids)
    exp_ids = organizer.load(files.virtual_experiment_probe_mesh_ids)

    for expFile in files.virtual_experiment_data_files:
        expData = organizer.load_red_file(expFile)
        probeData = expData[par.optimization_variable][exp_ids]

        opt_data = probeData[opt_probe_pnts]

        A = np.array([m.ndata[exp_name][opt_mesh_ids] for m in modes]).T


    for tankPresure, exp_data in organizer.load_selected_experiment_data():
        for exp_name in par.experiment_variable_mapping:
            exp_var_data = exp_data[exp_name]
            A = np.array([m.ndata[exp_name][opt_mesh_ids] for m in modes]).T
            coeffs = np.linalg.lstsq(A, exp_var_data[opt_probe_pnts])[0]
            organizer.save_coeff(tankPresure, exp_name, coeffs)


if __name__ == "__main__":
    from settings import *
    perform(dirs, files, par, organizer)