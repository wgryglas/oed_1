__author__ = 'wgryglas'

import collect_converged
import gather_all_sol_in_one_dir
import add_out_mach_to_data_filenames
import settings


modules = [add_out_mach_to_data_filenames] #collect_converged, gather_all_sol_in_one_dir, add_out_mach_to_data_filenames

for module in modules:
    module.perform(settings.dirs, settings.files, settings.par)

