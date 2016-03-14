import matplotlib.pyplot as plt
import numpy as np
from wg.tools.function import *
import modred as md

__author__ = 'wgryglas'


points = fixargument(plt.plot, marker='o', linestyle='', markersize=3)
thickLines = fixargument(plt.plot, marker='', linewidth=3)
plt.show = join(lambda : plt.get_current_fig_manager().full_screen_toggle(), plt.show)


plotnow = join(plt.figure, points, lambda title: plt.title(title), plt.show, argsFor={1: 'all', 2: 'title'})


#np.random.seed(1)

@assert_arg_iterable(args=0, kwargs='X')
def fun(X, a, b, c):
    res = []
    for x in X:
        res.append(a*x**2 + np.sin(b*x) + c)
        #res.append(b*np.sin(a*x**2) + c*x)
    return np.array(res)

# ==================================================================================================================== #
# Setup
# ==================================================================================================================== #
npoints = 50
n_mods = 3
max_measure = 6

# ==================================================================================================================== #
# Generate input data
# ==================================================================================================================== #
x = np.linspace(0, 3, npoints)
Y = np.zeros((npoints, 27))

count = 0
for i in range(1, 3):
    for j in range(1, 3):
        for k in range(1, 3):
            print i, j, k
            Y[:, count] = fun(x, i, j, k)
            count += 1

# ==================================================================================================================== #
# plot input data
# ==================================================================================================================== #
# f, ax = plt.subplots(1, 2)
# for i in range(count):
#     ax[0].plot(x, Y[:, i], 'o-', markersize=3)
#
# ax[0].grid(True)
# ax[1].grid(True)

# ==================================================================================================================== #
# Compute modes
# ==================================================================================================================== #
V, E = md.compute_POD_matrices_snaps_method(Y, range(n_mods))

# ==================================================================================================================== #
# plot modes
# ==================================================================================================================== #
# for i in range(n_mods):
#     ax[1].plot(x, V[:, i], 'o')
#
# ax[0].set_title('Rozwiazania dla kombinacji parametrow a,b,c=1,2')
# ax[1].set_title('Wektory bazowe(POD)')
#
# plt.show()

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
coeffs = [(1.2, 1.8, 1.5)]#(1, 2, 2), (2, 1, 1)]

exact = [fun(x, *c) for c in coeffs]
eMax = []
eMin = []
for e in exact:
    eMax.append(max(e))
    eMin.append(min(e))

eMax = max(eMax)
eMin = min(eMin)

rand = fixargument(np.random.normal, loc=0, scale=0.05*(eMax-eMin), size=1)

allMeasure = np.copy(exact)
for i in range(len(coeffs)):
    for j in range(npoints):
        allMeasure[i, j] = allMeasure[i, j] + rand()

# ==================================================================================================================== #
# Reconstruct and plot optimal
# ==================================================================================================================== #
# A = V[opt_pnts, :]
# dataX = x[opt_pnts]
#
# f, ax = plt.subplots(1, 2)
#
# legends=[]
# for e in exact:
#     legends.append(ax[1].plot(x, e, 'b--')[0])
#
# ax[1].vlines(dataX, eMin, eMax, linestyles='--', colors='gray')
#
# legend_str = ['Exact', 'Measurement', 'Reconstruction']
#
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
# ax[1].legend(legends, np.repeat(legend_str, len(legends)/3))

# ==================================================================================================================== #
# Reconstruct nad plot equal distribution
# ==================================================================================================================== #
# legends=[]
# for e in exact:
#     legends.append(ax[0].plot(x, e, 'b--')[0])
#
# n_measure = len(opt_pnts)
# di = int(npoints/(n_measure+1))
# dataId = [di*(i+1) for i in range(n_measure)]
# dataX = x[dataId]
#
# ax[0].vlines(dataX, eMin, eMax, linestyles='--', colors='gray')
#
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
# ax[0].legend(legends, np.repeat(legend_str, len(legends)/3))
#
# plt.show()



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













