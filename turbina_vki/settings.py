# coding=utf8
import sys
import os
import numpy as np
reload(sys)
sys.setdefaultencoding('UTF8')

__author__ = 'wgryglas'
"""
Settings file for VKI turbine cases
"""

# ------------------ DIRECTORIES ----------------------------------------
class __dirs__:
    """
    Class providing paths to certain directories required during
    optimization and data preparation.

    Changing root member variable automatically changes paths defined
    in 'dirs' and 'files' global variable.
    """
    def __init__(self, base): self.root = base

    @property
    def all_data(self): return self.root + os.sep + 'all'

    @property
    def converged_data(self): return self.root + os.sep + 'converged'

    @property
    def experiment_root(self): return self.root + os.sep + "experimental"

    @property
    def experiment(self): return self.experiment_root + os.sep + "Loptaki_23-03"

dirs = __dirs__('/home/wgryglas/AVIO/data/tunel_vki2')


# ------------------ FILES ----------------------------------------
class __files__:
    """
    Class providing paths to certain files required during
    optimization and data preparation
    """
    @property
    def boundary_coords(self): return dirs.root + os.sep + 'boundaryNodesCoordinates.dat'

    @property
    def probes_coords(self): return dirs.root + os.sep + "probe_points"

    @property
    def geom_get(self): return dirs.root + os.sep + "geom/kaskada_ns.get"

    @property
    def experiment_datas(self): return dirs.expriment + os.sep + np.array(["press0-675_pressp_1_pressp_2_machp1_machp2.npy", "press0-550_pressp_1_pressp_2_machp1_machp2.npy", "press0-225_pressp_1_pressp_2_machp1_machp2.npy"])

    @property
    def boundary_source(self): return dirs.root + os.sep + 'aoa-90.0ma1/_sol_surf.dat'

    @property
    def experiment_spreadsheet(self): return dirs.experiment_root + os.sep + "opracowanie_vki_ls-59.ods"


files = __files__()


# ------------------ PARAMETERS ----------------------------------------

class __parameters__:
    @property
    def kappa(self): return 1.4

    @property
    def experiment_sheet_name(self): return u'Łopatki 23.03'

    @property
    def red_output_file_name(self): return "_sol.dat"

    @property
    def experiment_output_suffix(self): return '_pressp_1_pressp_2_machp1_machp2'

par = __parameters__()


