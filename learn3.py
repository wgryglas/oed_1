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

npoints = 50
n_mods = 3
max_measure = 6

x = np.linspace(0, 3, npoints)
Y = np.zeros((npoints, 27))

count = 0
for i in range(1, 3):
    for j in range(1, 3):
        for k in range(1, 3):
            print i, j, k
            Y[:, count] = fun(x, i, j, k)
            count += 1

f, ax = plt.subplots(2, 2)

for i in range(count):
    ax[0, 0].plot(x, Y[:, i], '-', markersize=3)


V, E = md.compute_POD_matrices_snaps_method(Y, range(n_mods))

for i in range(n_mods):
    ax[0, 1].plot(x, V[:, i], '-')


# Optimization:
opt_pnts = [int(npoints/4), int(3*npoints/4)]

increase = list()

minimal = (0, float('inf'))
while minimal[0] != -1 and len(opt_pnts) < max_measure:
    minimal = (-1, minimal[1])
    locmin = float('inf')
    for j in range(npoints):
        if j in opt_pnts:
            continue

        tmp_pnts = [v for v in opt_pnts]
        tmp_pnts.append(j)
        W = V[tmp_pnts, :]
        M = W.T*W #W*W.T
        try:
            # tr = np.ravel(np.linalg.inv(M).trace())
            D = np.linalg.eigvalsh(M)
            tr = np.ravel(np.sum(1./D))
        except np.linalg.LinAlgError:
            continue

        if tr < 0 :
            continue

        #tr = np.abs(tr)

        if tr < locmin:
            locmin = tr

        if tr < minimal[1]:
            minimal = (j, tr)
        else:
            continue

    if minimal[0] != -1:
        print 'adding', minimal[0], 'with tr', minimal[1]
        opt_pnts.append(minimal[0])
    else:
        print 'not added, localmin=', locmin


# Generate measurement:
coeffs = [(1, 2, 2)]#, (2, 1, 1)]

exact = [fun(x, *c) for c in coeffs]
eMax = []
eMin = []
for e in exact:
    eMax.append(max(e))
    eMin.append(min(e))

eMax = max(eMax)
eMin = min(eMin)

rand = fixargument(np.random.normal, loc=0, scale=0.05*eMax, size=1)

allMeasure = np.copy(exact)
for i in range(len(coeffs)):
    for j in range(npoints):
        allMeasure[i,j] = allMeasure[i,j] + rand()

for e in exact:
    ax[1, 0].plot(x, e, 'b--')
    ax[1, 1].plot(x, e, 'b--')


# Reconstruct optimal
A = V[opt_pnts, :]
dataX = x[opt_pnts]

ax[1, 1].vlines(dataX, eMin, eMax, linestyles='--', colors='gray')

for measure in allMeasure:
    dataY = measure[opt_pnts]
    T = np.linalg.lstsq(A, dataY)[0]
    T = np.matrix([T])
    Ynew = V*T.T
    Ynew = np.ravel(Ynew)
    ax[1, 1].plot(dataX, dataY, 'ro')
    ax[1, 1].plot(x, Ynew, 'r-')


# # Reconstruct equal distribution
n_measure = len(opt_pnts)
di = int(npoints/(n_measure+1))
dataId = [di*(i+1) for i in range(n_measure)]
dataX = x[dataId]

ax[1,0].vlines(dataX, eMin, eMax, linestyles='--', colors='gray')

for measure in allMeasure:
    dataY = measure[dataId]
    A = V[dataId, :]
    T = np.linalg.lstsq(A, dataY)[0]
    T = np.matrix([T])
    Ynew = V*T.T
    Ynew = np.ravel(Ynew)
    ax[1, 0].plot(dataX, dataY, 'ro')
    ax[1, 0].plot(x, Ynew, 'r-')

plt.show()


# plt.figure()
# plt.plot([i for i in range(len(increase))], increase)
# plt.show()
