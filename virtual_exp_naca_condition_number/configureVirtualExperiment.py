
from add_pressure_mach_to_data import renameAddPressureMach
from bearded_octo_wookie.RED.RedReader import RedTecplotFile

from settings import par, dirs

f = dirs.experiment_dir + "/aoa2_2_ma_0_7.dat"
fout = dirs.experiment_dir + "/aoa2_2_ma_0_7-allvar.dat"
f_data = RedTecplotFile(f, useCache=True)


renameAddPressureMach(f_data, par.kappa, par.p_inf)

f_data.dumpFile(fout)