# coding=utf8
import sys
import os
import numpy as np
from wg.tools.system import ask_to_clean_dir
reload(sys)
sys.setdefaultencoding('UTF8')

__author__ = 'wgryglas'
"""
Settings file for VKI turbine cases
"""

# ------------------ PARAMETERS ----------------------------------------

class __parameters__:
    def __init__(self):
        self.num_modes = 5

        self.num_measure_pnts = 2 * self.num_modes + 1

        self.max_data_files = 200

        self.useCache = True

        self.optimization_variable = "mach_iso"

        self.experiment_variable_mapping = {"p": 0, "mach_iso": 2}

        self.plot = True#False

        self.save_plot_data = False#True


    @property
    def kappa(self): return 1.4

    @property
    def experiment_sheet_name(self): return u'Łopatki 23.03'

    @property
    def red_output_file_name(self): return "_sol.dat"

    @property
    def experiment_output_suffix(self): return '_pressp_1_pressp_2_machp1_machp2'

    @staticmethod
    def filter_calc_out_mach(mach):
        return 0.7 < mach #mach < 1.

    @staticmethod
    def filter_calc_angle_of_attack(angle):
        return True # 90 < angle < 91.5

    @staticmethod
    def filter_experiment_pressure_in_tank(pressure):
        return True #pressure >= 0.4


par = __parameters__()


# ------------------ DIRECTORIES ----------------------------------------
class __dirs__:
    """
    Class providing paths to certain directories required during
    optimization and data preparation.

    Changing root member variable automatically changes paths defined
    in 'dirs' and 'files' global variable.
    """
    def __init__(self, base):
        self.root = base

    @property
    def all_data(self): return self.root + os.sep + 'all'

    @property
    def converged_data(self): return self.root + os.sep + 'converged'

    @property
    def experiment_root(self): return self.root + os.sep + "experimental"

    @property
    def experiment(self): return self.experiment_root + os.sep + "Loptaki_23-03"

    @property
    def modes(self): return self.root + os.sep + "modes"

    @property
    def reconstructs(self): return self.root + os.sep + "reconstructs2"

    @property
    def results(self): return "set_"+str(1)

    @property
    def plot_data_save_dir(self): return '/home/wgryglas/Documents/studia doktoranckie/seminaria/za2017/figures/vki/sub_and_supersonic'


dirs = __dirs__('/home/wgryglas/AVIO/data/tunel_vki2b_v2')


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
    def probes_base_coords(self): return dirs.root + os.sep + "probe_base_points"

    @property
    def geom_get(self): return dirs.root + os.sep + "geom/kaskada_ns.get"

    @property
    def experiment_datas(self):
        from glob import glob1
        return glob1(dirs.experiment, "press*")

    @property
    def boundary_source(self): return dirs.root + os.sep + 'aoa-90.0ma1/_sol_surf.dat'

    @property
    def experiment_spreadsheet(self): return dirs.experiment_root + os.sep + "opracowanie_vki_ls-59.ods"

    @property
    def red_converged_results(self):
        from glob import glob1
        return [dirs.all_data + os.sep + f for f in glob1(dirs.all_data, "*.dat")]

    def modes(self, nmodes):
        """
        :return: list of tuples (mode index, mode file path)
        """
        return [(i, dirs.modes+os.sep+"mode_"+str(i)+".dat") for i in range(nmodes)]

    @property
    def modes_stored(self):
        from glob import glob1
        return [(i, dirs.modes+os.sep+"mode_"+str(i)+".dat") for i in range(len(glob1(dirs.modes, "mode_*.dat")))]

    @property
    def modes_coeffs_stored(self):
        from glob import glob1
        return [dirs.modes + os.sep + f for f in glob1(dirs.modes, "coeff*.npy")]

    def modes_coeff_for(self, pressure_in_tank, exp_data_name):
        return dirs.modes + os.sep + "coeff_"+(str(pressure_in_tank).replace(".", "-")+"_"+exp_data_name+".npy")

    def reconstructs_for(self, pressure_in_tank, opt_variable, num_modes):
        return dirs.reconstructs + os.sep + opt_variable + "_" + str(num_modes) + \
               "_"+(str(pressure_in_tank).replace(".", "-"))+".dat"

    def reconstructs_stored(self):
        from glob import glob1
        return [dirs.reconstructs + os.sep + f for f in glob1(dirs.reconstructs, "*.dat")]

    @property
    def eigen_values(self):
        return dirs.modes + os.sep + "eigen_values.npy"

    @property
    def optimized_mesh_ids(self): return dirs.modes + os.sep + "optimized_mesh_ids.npy"

    @property
    def optimized_probe_ids(self): return dirs.modes + os.sep + "optimized_probe_ids.npy"

    @property
    def inistial_mesh_ids(self): return dirs.modes + os.sep + "initial_mesh_ids.npy"

    @property
    def inistial_probe_ids(self): return dirs.modes + os.sep + "initial_probe_ids.npy"

    @property
    def pressure_in_tank_vs_outflow_mach(self): return dirs.experiment_root + os.sep + "pressuer_in_tank_vs_out_mach.txt"

    @property
    def probe_mesh_ids(self): return dirs.experiment + os.sep + "probe_mesh_ids.npy"


    # data for plotting:
    @property
    def plot_curve_xy(self): return dirs.plot_data_save_dir + os.sep + "curve_xy.csv"

    @property
    def plot_probe_pnts(self): return dirs.plot_data_save_dir + os.sep + "probe_pnts.csv"

    @property
    def plot_opt_pnts(self): return dirs.plot_data_save_dir + os.sep + "opt_pnts.csv"

    @property
    def plot_initial_pnts(self): return dirs.plot_data_save_dir + os.sep + "initial_pnts.csv"

    def plot_boundary_reconstruct(self, name): return dirs.plot_data_save_dir + os.sep + name + ".csv"





files = __files__()


# ------------------ DATA ORGANIZER ----------------------------------------
class __data_organizer__:
    def __init__(self):
        self.__red_cache__ = None
        self.__get_cache__ = None

    @staticmethod
    def read_angle_outmach(fname):
        try:
            s = fname.split("-")
            s2 = s[2].split("ma")[0]
            s = ".".join([s[1],s2])
            ang = float(s)
        except:
            ang = float(".".join(fname.split("-")[1]))
        return float(fname.split("outmach")[-1].split(".")[0].replace("-", ".")), ang

    @staticmethod
    def read_pressure_in_tank(fname):
        # import numpy as np
        # data = np.load(fname)
        # res = {}
        # for d in data:
        #     res[d[0]] = d[1]
        # return res
        pass

    @staticmethod
    def read_mode_id(fname):
        return int(fname.split("_")[1].split(".")[0])

    @staticmethod
    def select_mach_angle(fname):
        try:
            mach, angle = __data_organizer__.read_angle_outmach(os.path.basename(fname))
            return __parameters__.filter_calc_out_mach(mach) and __parameters__.filter_calc_angle_of_attack(angle)
        except:
            return False

    def load_selected_red_data(self):
        return self.load_red_data(__data_organizer__.select_mach_angle)

    def load_red_data(self, filterFun=None):
        from bearded_octo_wookie.RED.RedReader import RedTecplotFile

        if filterFun:
            prov = filter(filterFun, files.red_converged_results)
        else:
            prov = files.red_converged_results

        count = 0
        for f in prov:
            if count > par.max_data_files:
                raise StopIteration()

            yield RedTecplotFile(f, useCache=par.useCache)
            count += 1

    def load_modes(self):
        """
        Generator providing modes as tecplot loaded file
        :return: generator
        """
        from bearded_octo_wookie.RED.RedReader import RedTecplotFile
        for i, fpath in files.modes_stored:
            yield RedTecplotFile(fpath, useCache=par.useCache)

    def load_mode(self, modeNumber):
        from bearded_octo_wookie.RED.RedReader import RedTecplotFile
        for i, fpath in files.modes_stored:
            mid = __data_organizer__.read_mode_id(os.path.basename(fpath))
            if mid == modeNumber:
                return RedTecplotFile(fpath, useCache=par.useCache)
        return None

    def clear_cache(self, direcotry):
        from glob import glob1
        toRemove = [direcotry + os.sep + name.split(".")[0]+"cache.npy" for name in glob1(direcotry)]
        map(os.remove, toRemove)

    def clear_computed_data(self):
        from wg.tools.system import clear_dir
        clear_dir(dirs.modes)
        clear_dir(dirs.reconstructs)

    def clear_all_cache(self):
        self.clear_cache(dirs.all_data)
        self.clear_cache(dirs.modes)
        self.clear_cache(dirs.reconstructs)

    def load_mesh_boundary_coordinates(self):
        from RedBoundaryReader import readBoundaryNodes
        return readBoundaryNodes(files.boundary_coords)

    def get_mesh_boundary_ids(self):
        from scipy.spatial import cKDTree
        bX, bY = self.load_mesh_boundary_coordinates()
        for idata in self.load_red_data():
            return cKDTree(np.array([idata.ndata.X, idata.ndata.Y]).T).query(np.array([bX, bY]).T)[1]

    def load_experiment_points(self):
        eXY = np.load(dirs.experiment+os.sep+"points.npy")
        return eXY[:,0], eXY[:,1]

    def load_experiment_data(self, filterFun=None):
        """
        :return: list of tuples (pressure, 2 col data with pressure and mach)
        """
        for fname in files.experiment_datas:
            pressure = float(fname.split("_")[0][5:].replace("-", "."))

            if filterFun and not filterFun(pressure):
                continue

            data = np.load(dirs.experiment+os.sep+fname)
            exp_data = {}
            for n in par.experiment_variable_mapping:
                exp_data[n] = data[:, par.experiment_variable_mapping[n]]
            yield pressure, exp_data

    def load_selected_experiment_data(self):
        return self.load_experiment_data(__parameters__.filter_experiment_pressure_in_tank)

    def transform_to_experiment_csys(self, x, y):
        alpha = np.radians(-26.7)
        y = - y
        dx = 0.239928884974
        dy = 0.4564873
        bXNew = x*np.cos(alpha)-y*np.sin(alpha) + dx
        bYNew = x*np.sin(alpha)+y*np.cos(alpha) + dy
        return bXNew, bYNew

    def tranform_to_mesh_csys(self, x, y):
        alpha = np.radians(-26.7)
        dx = 0.239928884974
        dy = 0.4564873

        c = np.cos(alpha)
        s = np.sin(alpha)

        xNew = x*c + y*s - dx*c - dy*s
        yNew = x*s - y*c - dx*s + dy*c
        return xNew, yNew

    def save_coeff(self, pressure_in_tank, exp_data_name, coeffs):
        import numpy as np
        f = files.modes_coeff_for(pressure_in_tank, exp_data_name)
        np.save(f, coeffs)

    def load_mesh(self):
        """
        Loads mesh in cKDTree data structure.
        :return: cKDTree object constructed from mesh coordinates
        """
        from scipy.spatial import cKDTree
        inputData = self.load_selected_red_data().next()
        return cKDTree(inputData.data[:,:2])

    def load_coeffs(self):
        """
        :return: list of tuples (pressure_in_tank, variable name, array of coeffs)
        """
        import re

        for f in files.modes_coeffs_stored:
            data = np.load(f)
            groups = re.findall(r'coeff_([0-9]+-[0-9]*)_([0-9,aA-zZ]+).npy', os.path.basename(f))[0]
            press = float(groups[0].replace("-","."))
            varname = groups[1]
            yield press, varname, data

    def load_get_geometry(self):
        from bearded_octo_wookie.RED.RedReader import GetReader
        return GetReader(files.geom_get)

    def get_template_tecplot_file(self):
        if not self.__red_cache__:
            self.__red_cache__ = self.load_selected_red_data().next()
        return self.__red_cache__.copyForWriting()

    @property
    def number_of_modes_stored(self):
        return len(files.modes_stored)

    def sort_points(self, X, Y, curves=None):
        import numpy as np

        if not curves:
            if not self.__get_cache__:
                self.__get_cache__ = self.load_get_geometry()
            curves = self.__get_cache__

        XY = [[x, y] for x, y in zip(X, Y)]

        sort = np.argsort(curves.getParamList(XY, 1))
        return sort

    def save(self, path, data):
        import numpy as np
        np.save(path, data)


    def save_plot_data(self, path, data_dictionary):
        import pandas
        from wg.tools.system import savefile

        savefile(pandas.DataFrame(data_dictionary).to_csv, path, index=False)


    def load(self, path):
        import numpy as np
        return np.load(path)

    def get_reconstruct_info(self, fname):
        import re
        print fname
        info = re.findall(r'([0-9,aA-zZ]+[_[0-9,aA-zZ]*]*)_([0-9]+)_([0-9]+[-[0-9]+]*).dat', fname)[0]
        return info[0], int(info[1]), float(info[2].replace("-", "."))

    def load_reconstruction(self, filterFun=None):
        """
        :param filterFun:
        :return: tuple: (opt var, num modes, pressure in tank, tecplotfile)
        """
        import os
        from bearded_octo_wookie.RED.RedReader import RedTecplotFile

        for fpath in files.reconstructs_stored():
            fname = os.path.basename(fpath)
            infos = self.get_reconstruct_info(fname)
            if filterFun and not filterFun(*infos):
                continue
            yield infos[0], infos[1], infos[2], RedTecplotFile(fpath, useCache=par.useCache)

    def load_pressure_to_mach_mapping(self):
        import numpy as np
        return {p[0]:p[1] for p in np.loadtxt(files.pressure_in_tank_vs_outflow_mach)}


organizer = __data_organizer__()


# ------------------ ADDITIONAL TOOLS ---------------------------------------
def run_modules(*modules):
    """
    Function runs list of modules containing "perform" function or if provided
    object in list is the function itself then it runs this function.
    Each  callable object is called with 4 arguments: dirs, files, par, organizer
    :param module_list:
    :return: none
    """
    # decide if user provided list or var. num. of argunments
    if len(modules) == 1:
        module_list = modules[0]
    else:
        module_list = modules

    for module in module_list:
        if hasattr(module, "__call__"):
            print "--------------------------------------------"
            print "running function "+module.__name__ + " ..."
            print "--------------------------------------------"
            module(dirs, files, par, organizer)
        else:
            print "--------------------------------------------"
            print " running module "+module.__name__ + " ... "
            print "--------------------------------------------"
            module.perform(dirs, files, par, organizer)

