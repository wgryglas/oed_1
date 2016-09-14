

def perform(dirs, files, par, organizer):
    from bearded_octo_wookie.RED import POD
    from wg.tools.system import clear_dir

    idata = [d for d in organizer.load_selected_red_data()]

    nmodes = par.num_modes
    if len(idata) < nmodes:
        nmodes = len(idata)

    podmodes = POD.Modes(idata, num_modes=nmodes, appendAsVariables=False)

    clear_dir(dirs.modes)

    for mode, fpath in files.modes(nmodes):
        podmodes.dumpMode(fpath, mode)

    organizer.save(files.eigen_values, podmodes.eig_vals)

    return podmodes

if __name__ == "__main__":
    from settings import *
    perform(dirs, files, par, organizer)
