import glob, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

from bearded_octo_wookie.RED.RedReader import RedTecplotFile
from RedBoundaryReader import readBoundaryNodes
from tools import fit_experiment_porbe_locations_to_mesh_boundary
import matplotlib.animation as animation
from settings import *

__author__ = 'wgryglas'


kappa = 1.4

data_files = glob.glob1(dirs.all_data, "*.dat")

bX, bY = readBoundaryNodes(files.boundary_coords)
coords = RedTecplotFile(dirs.all_data+os.sep+data_files[0], useCache=True).data[:, :2]

mesh = cKDTree(coords)

bids = mesh.query(np.array([bX, bY]).T)[1]

p_inf = 1.


def read_and_compute(fname):
    f_data = RedTecplotFile(dirs.all_data+os.sep+fname, useCache=True).\
                 renameVariables(["X", "Y", "rho", "rhoE", "rhoU", "rhoV", "nut"]). \
                 computeVarialbe("p", lambda x, y, rho, e, u, v, nut: (kappa-1.) * (e - 0.5*(u**2 + v**2)/rho)).\
                 computeVarialbe("mach", lambda x, y, rho, e, u, v, nut, p: np.sqrt((u/rho)**2 + (v/rho)**2) / np.sqrt(kappa*p/rho)).\
                 computeVarialbe("mach_iso", lambda x, y, rho, e, u, v, nut, p, mach: np.sqrt(2./(1.4-1.)*((p_inf/p)**(0.4/1.4)-1.)))

    return f_data

files = glob.glob1(dirs.all_data, "*.dat")


def get_data():
    for fname in files:
        f_data = read_and_compute(fname)
        yield f_data.ndata


fig = plt.figure()
line, = plt.plot(mesh.data[bids, 0], mesh.data[bids, 0], '.')

eXY = np.load(dirs.expriment+os.sep+"points.npy")
pX, pY = fit_experiment_porbe_locations_to_mesh_boundary(eXY[:, 0], eXY[:, 1], bX, bY)
exp_data = np.ravel(np.load(dirs.expriment+os.sep+files.experiment_datas[0])[:, 2])
line2, = plt.plot(pX, exp_data/max(exp_data), '.')


# def distance(data1, data2):
#     return np.linalg.norm([min(data1 - d) for d in data2])
#
# min_dist = np.inf
# fid = -1
# for i, data in enumerate(get_data()):
#     d = distance(data.mach_iso[bids], exp_data)
#     if d < min_dist:
#        fid = i
#        min_dist = d
# line.set_ydata(read_and_compute(files[fid]).ndata.mach_iso[bids])
# plt.title(files[fid])
# plt.show()



# fname = "aoa-91-5ma0-75outmach0-737.dat"
# line.set_ydata(read_and_compute(fname).ndata.mach_iso[bids])
# plt.title(fname)
# plt.show()



gen = get_data()
def animate(t):
    data = gen.next()
    line.set_ydata((data.mach_iso[bids]))
    plt.gca().relim()
    plt.gca().autoscale_view()
    plt.title(files[t])
    return line,

anim = animation.FuncAnimation(fig, animate, frames=len(data_files), interval=500)

plt.show()


