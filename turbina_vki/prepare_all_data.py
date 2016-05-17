__author__ = 'wgryglas'

from settings import *
import collect_converged
import gather_all_sol_in_one_dir
import extract_experiment_data

# Set root directory to be processed
dirs.root = '/home/wgryglas/AVIO/data/tunel_vki'

# prepare all data
collect_converged.perform(dirs, files, par)
gather_all_sol_in_one_dir.perform(dirs, files, par)
extract_experiment_data.perform(dirs, files, par)

