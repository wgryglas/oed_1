

from settings import *

import compute_modes
import add_pressure_mach_to_data
import extract_boundary_nodes
import optimize_positions

# add_pressure_mach_to_data.perform(dirs, files, par, organizer)

# compute_modes.perform(dirs, files, par, organizer)

# extract_boundary_nodes.perform(dirs, files, par, organizer)

optimize_positions.perform(dirs, files, par, organizer)



