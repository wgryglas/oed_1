
from settings import dirs, files, par, run_modules, organizer
import collect_converged
import extract_boundary_nodes
import add_pressure_mach_to_data
import gather_all_sol_in_one_dir
import optimize_positions
import compute_modes
import plot_optimized_points
import extract_probe_ids_in_virtual_experiment_mesh
import extract_virtual_exp_boundary_nodes
import reconstruct_set_of_virtual_experiments

run_modules(reconstruct_set_of_virtual_experiments)