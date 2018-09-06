import matplotlib.pyplot as plt
import numpy as np
from wg.tools.function import *
from wg.tools.system import savefile
import modred as md

import wg.tools.system

from os import makedirs
from os.path import dirname
from os.path import exists


__author__ = 'wgryglas'

points = fixargument(plt.plot, marker='o', linestyle='', markersize=3)
thickLines = fixargument(plt.plot, marker='', linewidth=3)
#plt.show = join(lambda : plt.get_current_fig_manager().full_screen_toggle(), plt.show)
plt.figure = fixargument(plt.figure, figsize=(9, 4.5))

plotnow = join(plt.figure, points, lambda title: plt.title(title), plt.show, argsFor={1: 'all', 2: 'title'})

#def savefile(*args):
#    if not exists(dirname(args[0])):
#        makedirs(dirname(args[0]))
#    np.savetxt(*args)


#np.random.seed(1)

@assert_arg_iterable(args=0, kwargs='X')
def fun(X, a, b, c):
    res = []
    for x in X:
        res.append(a*x**2 + np.sin(b*x) +c)
        # res.append(b*np.sin(a*x**2) + c*x)
    return np.array(res)

# ==================================================================================================================== #
# Setup
# ==================================================================================================================== #
npoints = 50
n_mods = 5
max_measure =11#n_mods*2+1
n_div = 3

save_dir = "/home/wgryglas/Desktop/prezentacja_za/figures"#"/home/wgryglas/Documents/studia doktoranckie/seminaria/za2017"

# ==================================================================================================================== #
# Generate input data
# ==================================================================================================================== #
x = np.linspace(0, 3, npoints)
Y = np.zeros((npoints, n_div**3))

count = 0
for i in np.linspace(1, 3, n_div):
    for j in np.linspace(1, 3, n_div):
        for k in np.linspace(1, 3, n_div):
            print i, j, k
            Y[:, count] = fun(x, i, j, k)
            count += 1


data_save = np.zeros((npoints, n_div**3+1))
data_save[:, 0] = x
data_save[:, 1:] = Y
savefile(np.savetxt, save_dir+"/example1D/input.dat", data_save)

# ==================================================================================================================== #
# plot input data
# ==================================================================================================================== #
f = plt.figure()
ax = [f.add_subplot('121'), f.add_subplot('122')]
# f, ax = plt.subplots(1, 2)
for i in range(count):
    ax[0].plot(x, Y[:, i], 'o-', markersize=3)

for i in range(2):
    ax[i].set_xlabel("x")
    ax[i].set_ylabel("y")
    ax[i].grid(True)

# ==================================================================================================================== #
# Compute modes
# ==================================================================================================================== #
V, E = md.compute_POD_matrices_snaps_method(Y, range(n_mods))

# ==================================================================================================================== #
# plot modes
# ==================================================================================================================== #
for i in range(n_mods):
    ax[1].plot(x, V[:, i], 'o')

# ax[0].set_title('Rozwiazania dla kombinacji parametrow a,b,c=1,2')
# ax[1].set_title('Wektory bazowe(POD)')
#


#data_save=np.zeros((npoints, n_mods+1))
#data_save[:, 0] = x
#data_save[:, 1:] = V
#np.savetxt("/home/wgryglas/Desktop/prezentacja_za/example1D/basis.dat", data_save)

#plt.tight_layout()
#plt.savefig('/home/wgryglas/Desktop/prezentacja_esn/source_data.pdf', transparent=True)
plt.show()
# ==================================================================================================================== #
# Pnts location optimization:
# ==================================================================================================================== #
di = int(npoints/(max_measure+1))
opt_pnts = [di*(i+1) for i in range(max_measure)]

mintr = float('inf')
run = True
while run:
    for i in range(max_measure):
        tmp_opt = [v for v in opt_pnts]
        tmp_opt.pop(i)
        locmin = (-1, mintr)
        print 'moving ', i, 'point, tr=', mintr
        for j in range(npoints):
            if j not in tmp_opt:
                tmp_opt2 = [v for v in tmp_opt]
                tmp_opt2.append(j)
                W = V[tmp_opt2, :]
                M = W.T*W
                try:
                    D = np.linalg.eigvalsh(M)
                    tr = np.ravel(np.sum(1./D))
                except np.linalg.LinAlgError:
                    continue

                if tr < 0:
                    continue

                if tr < locmin[1]:
                    locmin = (j, tr)

        if locmin[0] != -1:
            opt_pnts[i] = locmin[0]
            mintr = locmin[1]
        else:
            run = False
            break

#opt_pnts  = [npoints/4, 3*npoints/4]
# minimal = (0, float('inf'))
# while minimal[0] != -1 and len(opt_pnts) < max_measure:
#     minimal = (-1, minimal[1])
#     locmin = float('inf')
#     for j in range(npoints):
#         if j in opt_pnts:
#             continue
#
#         tmp_pnts = [v for v in opt_pnts]
#         tmp_pnts.append(j)
#         W = V[tmp_pnts, :]
#         M = W.T*W #W*W.T
#         try:
#             # tr = np.ravel(np.linalg.inv(M).trace())
#             D = np.linalg.eigvalsh(M)
#             tr = np.ravel(np.sum(1./D))
#         except np.linalg.LinAlgError:
#             continue
#
#         if tr < 0 :
#             continue
#
#         #tr = np.abs(tr)
#
#         if tr < locmin:
#             locmin = tr
#
#         if tr < minimal[1]:
#             minimal = (j, tr)
#         else:
#             continue
#
#     if minimal[0] != -1:
#         print 'adding', minimal[0], 'with tr', minimal[1]
#         opt_pnts.append(minimal[0])
#     else:
#         print 'not added, localmin=', locmin

# ==================================================================================================================== #
# Generate measurement:
# ==================================================================================================================== #

nSamples = 100

coeffs = (1.2, 1.8, 1.5)#(1, 2, 2), (2, 1, 1)]

exact = fun(x, *coeffs)
eMax = max(exact)
eMin = min(exact)

rand = fixargument(np.random.normal, loc=0, scale=0.05*(eMax-eMin), size=1)

allMeasure = np.array([exact for i in range(nSamples)])
for i in range(nSamples):
    for j in range(npoints):
        allMeasure[i, j] = allMeasure[i, j] + rand()

# ==================================================================================================================== #
# Reconstruct and plot optimal
# ==================================================================================================================== #
A = V[opt_pnts, :]
dataX = x[opt_pnts]

f = plt.figure()
ax = [f.add_subplot('121'), f.add_subplot('122')]

# legends=[]
# for e in exact:
#     legends.append(ax[1].plot(x, e, 'b--')[0])

ax[1].vlines(dataX, eMin, eMax, linestyles='--', colors='gray')

recons = list()
for measure in allMeasure:
    dataY = measure[opt_pnts]
    T = np.linalg.lstsq(A, dataY)[0]
    T = np.matrix([T])
    Ynew = V*T.T
    Ynew = np.ravel(Ynew)
    recons.append(Ynew)

recons = np.array(recons)
y = np.mean(recons, axis=0)
stdderiv = np.std(recons, axis=0)
y1 = y - stdderiv
y2 = y + stdderiv
stdderiv_opt = stdderiv

legends=[]
legends.append(ax[1].plot(x, exact, 'b-', linewidth=1.5)[0])
legends.append(ax[1].plot(x, y, 'ko', marker='o', markersize=5, linewidth=1.5, alpha=1)[0])
legends.append(ax[1].plot(x, y1, color='k', linestyle='--', linewidth=0.5)[0])
ax[1].plot(x, y2, color='k', linestyle='--', linewidth=0.5)
ax[1].fill_between(x, y1, y2, where=y2 >= y1, facecolor='k', interpolate=True, alpha=0.2)
ax[1].set_xlim([-0.1, 3.1])
ax[1].set_ylim([1, 12])
ax[1].set_xlabel('x')
ax[1].set_ylabel('y')
# for measure in allMeasure:
#     dataY = measure[opt_pnts]
#     T = np.linalg.lstsq(A, dataY)[0]
#     T = np.matrix([T])
#     Ynew = V*T.T
#     Ynew = np.ravel(Ynew)
#     legends.append(ax[1].plot(dataX, dataY, 'ro')[0])
#     legends.append(ax[1].plot(x, Ynew, 'r-')[0])
#
# ax[1].set_title('Optimal positions')
ax[1].legend(legends, ['Exact', 'Reconstruction', 'Std. deviation'], loc=2)

savefile(np.savetxt, save_dir+"/alg1DCase-linearVsReconstruct/opt_pnts.dat", x[opt_pnts])

data_save = np.zeros((npoints, 6))
data_save[:, 0] = x
data_save[:, 1] = exact
data_save[:, 2] = y
data_save[:, 3] = stdderiv
data_save[:, 4] = y1
data_save[:, 5] = y2
savefile(np.savetxt, save_dir+"/alg1DCase-linearVsReconstruct/reconstruction.dat", data_save)

# ==================================================================================================================== #
# Reconstruct and plot equal distribution
# ==================================================================================================================== #
# legends=[]
# for e in exact:
#     legends.append(ax[0].plot(x, e, 'b--')[0])

n_measure = len(opt_pnts)
di = int(npoints/(n_measure-1))
dataId = [di*i for i in range(n_measure-1)]
dataId.append(npoints-1)

ax[0].vlines(x[dataId], eMin, eMax, linestyles='--', colors='gray')

from scipy.interpolate import interp1d

interpRecon = list()
for measure in allMeasure:
    inter = interp1d(x[dataId], measure[dataId], kind="linear")
    interpRecon.append(inter(x))

y = np.mean(interpRecon, axis=0)
stdderiv = np.std(interpRecon, axis=0)
y1 = y - stdderiv
y2 = y + stdderiv

legends=[]
legends.append(ax[0].plot(x, exact, 'b-', linewidth=2)[0])
legends.append(ax[0].plot(x, y, 'ko', markersize=5, linewidth=1)[0])
legends.append(ax[0].plot(x, y1, color='k', linestyle='--', linewidth=0.5)[0])
ax[0].plot(x, y2, color='k', linestyle='--', linewidth=0.5)
ax[0].fill_between(x, y1, y2, where=y2 >= y1, facecolor='k', interpolate=True, alpha=0.2)
ax[0].set_xlim([-0.1, 3.1])
ax[0].set_ylim([1, 12])
ax[0].set_xlabel('x')
ax[0].set_ylabel('y')
# for measure in allMeasure:
#     dataY = measure[dataId]
#     A = V[dataId, :]
#     T = np.linalg.lstsq(A, dataY)[0]
#     T = np.matrix([T])
#     Ynew = V*T.T
#     Ynew = np.ravel(Ynew)
#     legends.append(ax[0].plot(dataX, dataY, 'ro')[0])
#     legends.append(ax[0].plot(x, Ynew, 'r-')[0])
#
# ax[0].set_title('Equal distribution')
ax[0].legend(legends, ['Exact', 'L. interpolation', 'Std. deviation'], loc=2)


savefile(np.savetxt, save_dir+"/alg1DCase-linearVsReconstruct/eq_pnts.dat", x[dataId])

data_save = np.zeros((npoints, 6))
data_save[:, 0] = x
data_save[:, 1] = exact
data_save[:, 2] = y
data_save[:, 3] = stdderiv
data_save[:, 4] = y1
data_save[:, 5] = y2
savefile(np.savetxt, save_dir+"/alg1DCase-linearVsReconstruct/linear_interp.dat", data_save)

plt.tight_layout()
#plt.savefig('/home/wgryglas/Desktop/prezentacja_esn/recon_comparison.pdf', transparent=True)
plt.show()

plt.figure()
plt.legend([plt.plot(x, stdderiv, 'g-', marker='v')[0], plt.plot(x, stdderiv_opt, 'r-', marker='o')[0]], ['Interpolation', 'Reconstruction'])
plt.xlabel('x')
plt.ylabel('Std. deviation')
plt.grid()
#plt.savefig('/home/wgryglas/Desktop/prezentacja_esn/std_dev_comparison.pdf', transparent=True)
plt.show()
exit()


# ==================================================================================================================== #
# Interpolate modes
# ==================================================================================================================== #

import numpy.polynomial.chebyshev as cheb
import numpy.polynomial.legendre as leg
import numpy.polynomial.polynomial as pol

fitpol = pol.polyfit
evalpol = pol.polyval
#fitpol = cheb.chebfit
#evalpol = cheb.chebval
# fitpol = leg.legfit
# evalpol = leg.legval

x_fine = np.linspace(min(x), max(x), 2*npoints)

pols = []
# plt.figure()
for i in range(n_mods):
    pols.append( fitpol(x, V[:, i], V.shape[0]) )
    y_fine = np.ravel(cheb.chebval(x_fine, pols[i]))
    # plt.plot(x_fine, y_fine)
    # plt.plot(x, V[:, i], 'o')
# plt.show()


d = np.matrix(np.zeros((n_mods, n_mods)))
np.fill_diagonal(d, E)
dinv = np.matrix(np.zeros((n_mods, n_mods)))
np.fill_diagonal(dinv, 1./E)

U = np.matrix(Y).T * V * dinv

modes = np.zeros((2*npoints, n_mods))
for i in range(n_mods):
    modes[:, i] = np.ravel(evalpol(x_fine, pols[i]))

Res = modes*(d*U.T)

plt.figure()
for i in range(Y.shape[1]):
    plt.plot(x_fine, Res[:,i])
    plt.plot(x, Y[:,i], 'o')
plt._show()

# npoints = 20
# plt.figure()
# x = np.linspace(0, 2*np.pi, npoints)
# y = np.sin(x)
# x_fine = np.linspace(0, 2*np.pi, 2*npoints)
#
#
# pol = cheb.chebfit(x, y, V.shape[0])
# dPol = cheb.chebder(pol)
#
# y = np.ravel(cheb.chebval(x_fine, pol))
# dy = np.ravel(cheb.chebval(x_fine, dPol))
# dy_ex = np.cos(x_fine)
#
# plt.plot(x_fine, y)
# plt.plot(x_fine, dy)
# plt.plot(x_fine, dy_ex, 'o')
# plt.show()













