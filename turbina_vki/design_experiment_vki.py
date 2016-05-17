__author__ = 'wgryglas'

import glob, os, re
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

import modred as md
from bearded_octo_wookie.RED.RedReader import RedTecplotFile
from bearded_octo_wookie.RED import POD
from bearded_octo_wookie.RED.RedReader import GetReader
from RedBoundaryReader import readBoundaryNodes
from tools import fit_experiment_porbe_locations_to_mesh_boundary
from wg.tools.array import *
from wg.tools.function import join

from settings import *

plt.show = join(lambda : plt.get_current_fig_manager().full_screen_toggle(), plt.show)

# ------- SETUP ------------------------------------------------------------------------------------------------------ #
num_modes = 10
num_test_points = 2 * num_modes + 1
testVariable = 6 #<-("rho", "rhoE", "rhoU", "rhoV", "nut", "p", "mach_iso")
maxFiles = 100
exp_variable = 2 # <-- (pressure palisade 1, pressure palisade 2, mach palisade 1, mach palisade 2
p_inf = 1.
# ------- PREPARE DATA ----------------------------------------------------------------------------------------------- #
inputData = list()
for i, fname in enumerate(glob.glob1(dirs.all_data, "*.dat")):
    if i > maxFiles:
        break

    inputData.append(RedTecplotFile(dirs.all_data+os.sep+fname, useCache=True).
                     renameVariables(["X", "Y", "rho", "rhoE", "rhoU", "rhoV", "nut"]).
                     computeVarialbe("p", lambda x, y, rho, e, u, v, nut: (par.kappa-1.) * (e - 0.5*(u**2 + v**2)/rho)).
                     computeVarialbe("mach_iso", lambda x, y, rho, e, u, v, nut, p: np.sqrt(np.abs(2./(par.kappa-1.)*((p_inf/p)**((par.kappa-1)/par.kappa)-1.)))))
                     #computeVarialbe('p', '0.4*(rho-0.5*(rhoU**2+rhoV**2)/rho)'))


# read boundary nodes
mesh = cKDTree(inputData[0].data[:, :2])

bX, bY = readBoundaryNodes(files.boundary_coords)

boundary_ids = mesh.query(np.array([bX, bY]).T)[1]

# read geometry
curves = GetReader(files.geom_get)

# read probe points
#pX, pY = readBoundaryNodes(probes_coords_file)
eXY = np.load(dirs.expriment+os.sep+"points.npy")
pX, pY = fit_experiment_porbe_locations_to_mesh_boundary(eXY[:, 0], eXY[:, 1], bX, bY)


# find closest indices of boundary nodes to the probe coords
bProbeIds = cKDTree(np.array([bX, bY]).T).query(np.array([pX, pY]).T)[1]
probe_mesh_ids = boundary_ids[bProbeIds]

pX2 = mesh.data[probe_mesh_ids, 0]
pY2 = mesh.data[probe_mesh_ids, 1]

# find boundary and probe nodes parametrization (for setting appropriate order)

bParams = curves.getParamList(mesh.data[boundary_ids, :].tolist(), 1)
boundary_ids = boundary_ids[np.argsort(bParams)]

pParams = curves.getParamList(mesh.data[probe_mesh_ids, :].tolist(), 1)
probe_sort = np.argsort(pParams)
probe_mesh_ids = probe_mesh_ids[probe_sort]


# # Test ordering
# plt.figure()
# plt.plot(mesh.data[probeIds, 0], mesh.data[probeIds, 1], '.')
# for i, p in enumerate(mesh.data[probeIds, :]):
#     plt.annotate(str(i),p, p+[0.05, 0.05], arrowprops=dict(arrowstyle='->'))
# plt.show()


# # Plot initial points in ref. to all possible points:
# plt.figure()
# plt.plot(mesh.data[startIds, 0], mesh.data[startIds, 1], 'g.', markersize=20)
# plt.plot(bX, bY, 'g.', markersize=1)
# plt.plot(mesh.data[probeIds, 0], mesh.data[probeIds, 1], 'r.')
# plt.show()


def getSubFisherMatrix(modes, pntIds):
    X = np.matrix(modes[pntIds, :])
    return X.T*X

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


# Perform optimization

# Take initial positions:
start_mesh_ids = [i for i in spread(probe_mesh_ids, num_test_points)]
start_probe_ids = [i for i in spread_id(probe_mesh_ids, num_test_points)]

opt_pnts = np.copy(start_probe_ids)
n_pnts = len(start_probe_ids)

# Compute POD modes
podmodes = POD.Modes(inputData, num_modes=num_modes)

podmodes.plot_spectrum()

maxTest = 5
mintr = float('inf')

# Calc appropriate input data
varmodes = np.array([np.ravel(m[:, testVariable]) for m in podmodes.modes]).T

for loopId in range(maxTest):

    moved = False

    for i in range(n_pnts):
        mask = np.ones(n_pnts, dtype=bool)
        mask[i] = False
        tmp_opt = opt_pnts[mask]

        locmin = (-1, mintr)

        print 'checking', i, 'point with boundary position ', opt_pnts[i], 'and tr=', mintr

        traces = list()

        for j, _ in enumerate(probe_mesh_ids):
            if j in tmp_opt:
                continue

            tmp_opt2 = np.append(tmp_opt, j)

            testPos = probe_mesh_ids[tmp_opt2]
            criterion = computeACriterion(getSubFisherMatrix(varmodes, testPos))

            if criterion < 0:
                continue

            if criterion < locmin[1]:
                locmin = (j, criterion)

        if locmin[0] != -1:
            opt_pnts[i] = locmin[0]
            mintr = locmin[1]
            moved = True
            print 'moved to', locmin[0], 'tr=', locmin[1]

    print '---- END OF LOOP', loopId, '-----------------------------------------------------'

    if not moved:
        break


opt_mesh_ids = probe_mesh_ids[opt_pnts]



# Plot result of optimization
plt.figure()
plt.plot(mesh.data[opt_mesh_ids, 0], mesh.data[opt_mesh_ids, 1], 'r.', markersize=30)
plt.plot(mesh.data[start_mesh_ids, 0], mesh.data[start_mesh_ids, 1], 'g.', markersize=20)
plt.plot(mesh.data[probe_mesh_ids, 0], mesh.data[probe_mesh_ids, 1], 'b.')
plt.show()


#Reconstruct modes basing on numerical data:
# A = varmodes[opt_mesh_ids, :]
# in_file_id = 99
# coeffs = np.linalg.lstsq(A, inputData[in_file_id].data[:, testVariable][probe_mesh_ids][opt_pnts])[0]
# plt.figure()
# recvar = np.dot(varmodes, coeffs)
# plt.plot(mesh.data[boundary_ids,0], recvar[boundary_ids], '.')
# plt.plot(mesh.data[boundary_ids,0], inputData[in_file_id].data[:, testVariable][boundary_ids], '.')
# plt.show()


# Reconstruct fields
 # load experiment data and results in curve parameter order

reconstructs = list()
exp_datas = list()
for exp_data_file in files.experiment_datas:
    exp_data = np.ravel(np.load(dirs.expriment+os.sep+exp_data_file)[probe_sort, exp_variable])
    exp_datas.append(exp_data)
    A = varmodes[opt_mesh_ids, :]
    coeffs = np.linalg.lstsq(A, exp_data[opt_pnts])[0]
    reconstructs.append(np.dot(varmodes, coeffs))


def shift_color(col, factor):
    ncol = [c*(min([max([1.-factor, 0]), 1])) for c in col]
    ncol[3] = col[3]
    print ncol
    return ncol


# cmap = plt.get_cmap('jet_r')
for i, (rec, exp_data, fname) in enumerate(zip(reconstructs, exp_datas, files.experiment_datas)):
    plt.figure()
    value = re.search(r"[0-9]-[0-9]+", fname).group(0).replace('-', '.')
    #col = cmap(float(i)/len(exp_datas))
    col = [1., 0., 0., 1.]
    plt.plot(mesh.data[boundary_ids, 0], rec[boundary_ids], '.', color=col, markersize=5)
    plt.plot(pX[opt_pnts], exp_data[opt_pnts], '.', color="yellow", markersize=20)
    plt.plot(pX, exp_data, '.', color=shift_color(col, 0.25), markersize=15)
    plt.title("Pressure "+value)
    plt.show()
    # plt.plot(pX[opt_pnts], exp_data[opt_pnts], '.', color=col, markersize=10)
