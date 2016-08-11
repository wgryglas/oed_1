
def perform(dirs, files, par, organizer):

    from wg.tools.system import assert_dir_exists
    from wg.tools.system import clear_dir
    import numpy as np

    modes = [m for m in organizer.load_modes()]

    datas = {}

    import matplotlib.pyplot as plt
    bX, _ = organizer.load_mesh_boundary_coordinates()
    bids = organizer.get_mesh_boundary_ids()
    pX, _ = organizer.tranform_to_mesh_csys(*organizer.load_experiment_points())
    exp={press:data for (press, data) in organizer.load_selected_experiment_data()}

    for press, var, coeff in organizer.load_coeffs():
        print press, var

        varmodes = np.array([m.ndata[var] for m in modes]).T

        if press not in datas:
            datas[press] = organizer.get_template_tecplot_file()

        result = varmodes.dot(coeff)

        datas[press].appendData(var, result)

    clear_dir(dirs.reconstructs)

    for press in datas:
        datas[press].dumpFile(files.reconstructs_for(press, par.optimization_variable, par.num_modes))


if __name__ == "__main__":
    from settings import *
    perform(dirs, files, par, organizer)


