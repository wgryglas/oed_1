from scipy.stats.distributions import expon_gen

__author__ = 'wgryglas'

"""
Script checking import from experiment
"""

from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import glob
from bearded_octo_wookie.RED.POD import GetReader
from tools import fit_experiment_porbe_locations_to_mesh_boundary
from bearded_octo_wookie.RED.RedReader import RedTecplotFile
from RedBoundaryReader import readBoundaryNodes
from settings import *

data_file = dirs.all_data+os.sep+glob.glob1(dirs.all_data, "*.dat")[0]

data = RedTecplotFile(data_file, useCache=True)

mX = data.data[:, 0]
mY = data.data[:, 1]

bX, bY = readBoundaryNodes(files.boundary_coords)

mesh = cKDTree(data.data[:, :2])

boundary_mesh_ids = mesh.query(np.array([bX, bY]).T)[1]


eXY = np.load(dirs.experiment + os.sep+"points.npy")
eX, eY = fit_experiment_porbe_locations_to_mesh_boundary(eXY[:, 0], eXY[:, 1], bX, bY)
probe_boundary_ids = cKDTree(np.array([bX, bY]).T).query(np.array([eX, eY]).T)[1]

# PLOT probes over boundary
plt.figure()
plt.plot(eX, eY, 'r.', markersize=10)
plt.plot(mX[boundary_mesh_ids[probe_boundary_ids]], mY[boundary_mesh_ids[probe_boundary_ids]], 'g.', markersize=10)
plt.show()


# read exp. data
exp_data = np.load(files.experiment_datas[0])

# read geometry
curves = GetReader(files.geom_get)

probe_curve_param = curves.getParamList(np.array([bX[probe_boundary_ids], bY[probe_boundary_ids]]).T.tolist(), 1)

# PLOT data over parameter space
plt.figure()
# plt.plot(probe_curve_param, exp_data[:, 0], '.')
plt.plot(eX, exp_data[:, 0], '.')
# for i, (p, d) in enumerate(zip(probe_curve_param, exp_data[:, 0])):
#     plt.annotate(str(i), (p,d), (p+0.05, d))
plt.show()