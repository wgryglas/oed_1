

from settings import dirs, files, par, run_modules, organizer
import collect_converged
import extract_boundary_nodes
import add_pressure_mach_to_data
import gather_all_sol_in_one_dir
import optimize_positions

# run_modules([optimize_positions])


geom = organizer.load_get_geometry()
import numpy as np


xy = np.array( geom.getXYList(np.linspace(0, 1, 100), 1) )

import matplotlib.pyplot as plt

bX, bY = organizer.load_mesh_boundary_coordinates()

plt.plot(xy[:,0], xy[:,1], '.')
plt.plot(bX, bY, ".k", markersize=5)
plt.show()

