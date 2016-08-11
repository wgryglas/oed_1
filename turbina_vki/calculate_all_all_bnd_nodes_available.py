
import compute_modes
import optimize_all_positions
import fit_coefficent_to_experiment
import reconstruct_fields

import settings

# settings.run_modules(compute_modes, optimize_all_positions)
optimize_all_positions.perform(settings.dirs, settings.files, settings.par, settings.organizer)



