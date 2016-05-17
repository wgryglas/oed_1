__author__ = 'wgryglas'


from bearded_octo_wookie.RED.RedReader import RedTecplotFile
from bearded_octo_wookie.RED import POD
from bearded_octo_wookie.RED.RedReader import GetReader
from RedBoundaryReader import readBoundaryNodes
from tools import fit_experiment_porbe_locations_to_mesh_boundary
from wg.tools.array import *
from wg.tools.function import join

import numpy as np
import matplotlib.pyplot as plt
kappa = 1.4

data = RedTecplotFile("/home/wgryglas/AVIO/data/tunel_vki/aoa-85.5ma0.35/_sol_surf_trimmed.dat", useCache=True).\
                     renameVariables(["X", "Y", "VnX", "VnY",  "rho", "rhoE", "rhoU", "rhoV", "nut"]).\
                     computeVarialbe("p", lambda x, y, vnx, vny, rho, e, u, v, nut: (kappa-1.) * (e - 0.5*(u**2 + v**2)/rho)).\
                     computeVarialbe("mach", lambda x, y, vnx, vny, rho, e, u, v, nut, p: np.sqrt( (u/rho)**2 + (v/rho)**2) / np.sqrt(kappa*p/rho))



plt.figure()
plt.plot(data.ndata.X, data.ndata.mach)
plt.show()