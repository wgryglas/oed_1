import collect_converged
import gather_all_sol_in_one_dir
import add_out_mach_to_data_filenames
import extract_boundary_nodes
import extract_experiment_data
import settings
import add_pressure_mach_to_data


settings.run_modules( extract_experiment_data,
                      extract_boundary_nodes,
                      collect_converged,
                      gather_all_sol_in_one_dir,
                      add_pressure_mach_to_data,
                      add_out_mach_to_data_filenames,
                      lambda d, f, p, o: o.clear_all_cache)


