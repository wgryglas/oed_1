

from settings import *

import compute_modes
import add_pressure_mach_to_data
import extract_boundary_nodes
import optimize_positions
import plot_optimized_points
import reconstruct_set_of_virtual_experiments
import extract_probe_ids_in_virtual_experiment_mesh

# add_pressure_mach_to_data.perform(dirs, files, par, organizer)

# compute_modes.perform(dirs, files, par, organizer)

# extract_boundary_nodes.perform(dirs, files, par, organizer)

# optimize_positions.perform(dirs, files, par, organizer)


# plot_optimized_points.perform(dirs, files, par, organizer)

# run_modules([extract_probe_ids_in_virtual_experiment_mesh])

run_modules([reconstruct_set_of_virtual_experiments])#compute_modes, optimize_positions, plot_optimized_points, reconstruct_set_of_virtual_experiments
