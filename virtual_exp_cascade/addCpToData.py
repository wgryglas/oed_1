def addCpToData(redData, p_inf):
    data = redData.ndata

    refPoint = np.argmax(data.Y)

    ref_p_inf = data.p[refPoint]
    refDynamicPressure = (data.rhoU[refPoint]**2 + data.rhoV[refPoint]**2) / data.rho[refPoint] / 2

    redData.computeVarialbe("cp", lambda x, y, rho, e, u, v, nut, p, mach, mach_iso, cp: (p - ref_p_inf) / refDynamicPressure)

def perform(dirs, files, par, organizer):
    from bearded_octo_wookie.RED.RedReader import RedTecplotFile

    for f in files.red_converged_results:
            print "Processing file", f
            f_data = RedTecplotFile(f, useCache=False)
            addCpToData(f_data, par.p_inf)
            f_data.dumpFile(f)

    for f in files.virtual_experiment_data_files:
        print "Procssing v.e. data", f
        f_data = RedTecplotFile(f, useCache=False)
        addCpToData(f_data, par.p_inf)
        f_data.dumpFile(f)

    # Force program to load new, modified data in the new run
    organizer.clear_all_cache()


if __name__ == "__main__":
    from settings import *
    perform(dirs, files, par, organizer)
