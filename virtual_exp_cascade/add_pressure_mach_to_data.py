def renameAddPressureMach(redData, kappa, p_inf):
    import numpy as np
    powKappa = (kappa - 1.)/kappa
    redData.renameVariables(["X", "Y", "rho", "rhoE", "rhoU", "rhoV", "nut"])
    redData.computeVarialbe("p", lambda x, y, rho, e, u, v, nut: (kappa-1.) * (e - 0.5*(u**2 + v**2)/rho))
    redData.computeVarialbe("mach", lambda x, y, rho, e, u, v, nut, p: np.sqrt( (u/rho)**2 + (v/rho)**2) / np.sqrt(kappa*p/rho))
    redData.computeVarialbe("mach_iso", lambda x, y, rho, e, u, v, nut, p, mach: np.sqrt( np.abs( 2./(kappa -1. ) * ( np.abs(p_inf/p)**(powKappa) - 1. ) ) ) )

def perform(dirs, files, par, organizer):
    from bearded_octo_wookie.RED.RedReader import RedTecplotFile
    p_inf = 1.
    powKappa = (par.kappa - 1.)/par.kappa

    for f in files.red_converged_results:

            f_data = RedTecplotFile(f, useCache=False)

            if f_data.variables[0] == "X" and f_data.variables[-1]=="mach_iso":
                print "ommiting "+f
                continue

            print "processing "+f

            renameAddPressureMach(f_data, par.kappa, p_inf)

            f_data.dumpFile(f)




if __name__ == "__main__":
    from settings import *
    perform(dirs, files, par)