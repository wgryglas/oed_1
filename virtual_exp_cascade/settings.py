# coding=utf8
import sys
import os
import numpy as np
import itertools

from wg.tools.system import ask_to_clean_dir

reload(sys)
sys.setdefaultencoding('UTF8')

"""
Settings file for cascade
"""

def computeACriterion(M):
        try:
            D = np.linalg.eigvalsh(M)
            tr = np.ravel(np.sum(1./D))
            if tr < 0:
                tr = np.inf
        except Exception as e:
            print e
            tr = np.inf
        return tr

def computeConditionNumber(M):
    return np.linalg.cond(M)


# ------------------ PARAMETERS ----------------------------------------

class __parameters__:
    def __init__(self):
        self.num_modes = 10
        self.num_measure_pnts = 2 * self.num_modes + 1
        self.max_data_files = 1e4
        self.useCache = True
        self.optimization_variable = "cp" #"mach_iso" #"p"
        self.standard_deviation = 0.05
        self.virtual_exp_distortion_iteration = 200
        self.plot = False
        self.save_plot_data = True
        self.plot_in_parameter_space = False

        self.criterionFunction = computeACriterion
        # self.criterionFunction = computeConditionNumber

    @property
    def kappa(self): return 1.4

    @property
    def p_inf(self): return 1.

    @property
    def red_output_file_name(self): return "_sol.dat"



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
    def modes(self): return self.root + os.sep + "modes"

    @property
    def experiment_dir(self): return self.root + os.sep + "virtual_experiment"

    @property
    def reconstructs(self): return self.root + os.sep + "reconstructs"

    @property
    def results(self): return "set_"+str(1)

    @property
    def figures(self): return self.root + "/out/figures"

    @property
    def save_plot_dir(self): return "/home/wgryglas/Documents/studia doktoranckie/articles/COOPERNIK-Instrumentation/figures/cascade"
     #"/home/wgryglas/Desktop/figures_cascade_data"
     # #"/home/wgryglas/Documents/studia doktoranckie/articles/doe_pod_paper/figures/cascade/new_data"


dirs = __dirs__('/home/wgryglas/AVIO/data/cascade')


# ------------------ FILES ----------------------------------------
class __files__:
    """
    Class providing paths to certain files required during
    optimization and data preparation
    """
    @property
    def modes_test_data_file(self): return {r'\noindent Anagle of attack $0^0$\newline Outlet pressure $0.6$': dirs.all_data + os.sep + "aoa-90-0p0-6.dat",
                                            r'\noindent Anagle of attack $10^0$\newline Outlet pressure $0.8$': dirs.all_data + os.sep + "aoa-100-5p0-8.dat",
                                            r'\noindent Anagle of attack $25^0$\newline Outlet pressure $0.95$': dirs.all_data + os.sep + "aoa-125p0-95.dat"}

    @property
    def boundary_coords(self): return dirs.root + os.sep + 'boundary.dat'

    @property
    def probes_coords(self): return dirs.root + os.sep + "probe_points"

    @property
    def probes_base_coords(self): return dirs.root + os.sep + "probe_base_points"

    @property
    def geom_get(self): return dirs.root + os.sep + "geom/geometry.get"

    @property
    def boundary_source(self): return dirs.root + os.sep + 'aoa-100.5p0.6/_sol_surf.dat'

    @property
    def virtual_experiment_boundary_source(self): return dirs.experiment_dir + os.sep + '_out_surf.dat'

    @property
    def virtual_experiment_boundary_coords(self): return dirs.experiment_dir + os.sep + 'boundary.dat'

    @property
    def virtual_experiment_data_files(self): return [dirs.experiment_dir + os.sep + '_out-allvar.dat']

    @property
    def virtual_experiment_probe_mesh_ids(self): return dirs.experiment_dir + os.sep + 'probe_mesh_ids.npy'

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
    def probe_mesh_ids(self): return dirs.root + os.sep + "probe_mesh_ids.npy"


    @property
    def plot_profile(self): return dirs.save_plot_dir + os.sep + "profile.csv"

    @property
    def plot_init_opt_pnts(self): return dirs.save_plot_dir + os.sep + ("points_np%d_nm%d.csv" % (par.num_measure_pnts, par.num_modes)) # "points.csv"

    def plot_exact_result(self, data_name): return dirs.save_plot_dir + os.sep + data_name + os.sep + "exact.csv"

    def plot_reconstruction_result(self, data_name): return dirs.save_plot_dir + os.sep + data_name + os.sep + "reconstructionResult.csv"

    def plot_reconstruction_result_in_measurement(self, data_name): return dirs.save_plot_dir + os.sep + data_name + os.sep + "reconstructionResult-measurement_positions.csv"

    def plot_reconstruction_lowerCurve(self, data_name): return dirs.save_plot_dir + os.sep + data_name + os.sep + "reconstructionResult-lowerCurve.csv"

    def plot_reconstruction_upperCurve(self, data_name): return dirs.save_plot_dir + os.sep + data_name + os.sep + "reconstructionResult-upperCurve.csv"

    def plot_lininterpolation_result(self, data_name): return dirs.save_plot_dir + os.sep + data_name + os.sep + "linInterpResult.csv"

    def plot_lininterpolation_result_in_measurement(self, data_name): return dirs.save_plot_dir + os.sep + data_name + os.sep + "linInterpResult-measurement_positions.csv"

    def plot_lininterpolation_lowerCurve(self, data_name): return dirs.save_plot_dir + os.sep + data_name + os.sep + "linInterpResult-lowerCurve.csv"

    def plot_lininterpolation_upperCurve(self, data_name): return dirs.save_plot_dir + os.sep + data_name + os.sep + "linInterpResult-upperCurve.csv"



files = __files__()


# ------------------ DATA ORGANIZER ----------------------------------------
class __data_organizer__:
    def __init__(self):
        self.__red_cache__ = None
        self.__get_cache__ = None
        self.markers = itertools.cycle(['^', 'o', ' D', 's'])

    @staticmethod
    def read_mode_id(fname):
        return int(fname.split("_")[1].split(".")[0])

    def load_selected_red_data(self):
        return self.load_red_data()

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

    def load_red_file(self, path):
        from bearded_octo_wookie.RED.RedReader import RedTecplotFile
        return RedTecplotFile(path, useCache=par.useCache)

    def load_modes(self, numberOfModesToLoad = np.inf):
        """
        Generator providing modes as tecplot loaded file
        :return: generator
        """
        from bearded_octo_wookie.RED.RedReader import RedTecplotFile
        for i, fpath in files.modes_stored:
            yield RedTecplotFile(fpath, useCache=par.useCache)
            if i >= numberOfModesToLoad:
                return

    def load_mode(self, modeNumber):
        from bearded_octo_wookie.RED.RedReader import RedTecplotFile
        for i, fpath in files.modes_stored:
            mid = __data_organizer__.read_mode_id(os.path.basename(fpath))
            if mid == modeNumber:
                return RedTecplotFile(fpath, useCache=par.useCache)
        return None

    def clear_cache(self, directory):
        from glob import glob1
        toRemove = [directory + os.sep + name for name in glob1(directory, "*.cache.npz")]
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
        data = self.load_mesh_boundary_coordinates()
        return self.load_mesh().query(np.array(data).T)[1]

    def load_experiment_points(self):
        eXY = self.load_data_file(files.virtual_experiment_boundary_coords)
        return eXY[:, 0], eXY[:, 1]

    def load_virtual_exp_mesh(self):
        red = self.load_red_file(files.virtual_experiment_data_files[0])
        from scipy.spatial import cKDTree
        return cKDTree(red.data[:, :2])

    def load_virtual_exp_data(self):
        for f in files.virtual_experiment_data_files:
            yield self.load_red_file(f)

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

        if not self.__get_cache__:
            self.__get_cache__ = GetReader(files.geom_get)

        return self.__get_cache__

    def get_equaly_distributed_points_over_boundary(self, num_pnts):
        geom = self.load_get_geometry()
        return geom.getXYList(np.linspace(0, 1, num_pnts), 1)


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

    def load_data_file(self, path):
        import numpy as np
        return np.loadtxt(path)

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
        return {p[0]: p[1] for p in np.loadtxt(files.pressure_in_tank_vs_outflow_mach)}

    def get_probe_positions(self):
        return self.load_mesh_boundary_coordinates()

    def get_figure_setup(self):
        return {"figsize":(10, 7), "dpi": 360}

    def get_big_figure_setup(self):
        return {"figsize":(15, 8), "dpi": 360}

    @property
    def next_marker(self):
        return self.markers.next()


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
    if hasattr(modules, "__iter__"):
        module_list = modules
    else:
        module_list = [modules]

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

