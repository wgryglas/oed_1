import matplotlib.pyplot as plt
import numpy as np
from wg.tools.function import *
import modred as md

__author__ = 'wgryglas'


points = fixargument(plt.plot, marker='o', linestyle='', markersize=3)
thickLines = fixargument(plt.plot, marker='', linewidth=3)
plt.show = join(lambda : plt.get_current_fig_manager().full_screen_toggle(), plt.show)


plotnow = join(plt.figure, points, lambda title: plt.title(title), plt.show, argsFor={1: 'all', 2: 'title'})


@assert_arg_iterable(args=0, kwargs='X')
def fun(X, a, b, c):
    res = []
    for x in X:
        #res.append(a*x**2 + b*x + c)
        res.append(b*np.sin(a*x**2)+c*x)
    return np.array(res)

npoints = 21
n_mods = 3
n_mesure = 8

x = np.linspace(-1, 1, npoints)
Y = np.zeros((npoints, 27))

count = 0
for i in range(1, 3):
    for j in range(1, 3):
        for k in range(1, 3):
            print i, j, k
            Y[:, count] = fun(x, i, j, k)
            count += 1

f, ax = plt.subplots(1, 3)

for i in range(count):
    ax[0].plot(x, Y[:, i], '-o', markersize=3)


V, E = md.compute_POD_matrices_snaps_method(Y, range(n_mods))

for i in range(n_mods):
    ax[1].plot(x, V[:, i], '-')


coeffs = (1.2, 1.8, 1.5)
rand = fixargument(np.random.normal, loc=0, scale=0.05, size=1)

print n_mesure, npoints
di = int((npoints)/(n_mesure))
print di

dataId = [di*i for i in range(n_mesure)]
dataX = x[dataId]

dataY = fun(dataX, *coeffs)
for i in range(len(dataY)):
    dataY[i] = dataY[i] + rand()

ax[2].plot(dataX, dataY, 'ro')


exact = fun(x,*coeffs)
ax[2].plot(x, exact, 'b--')


A = V[dataId, :]

T = np.linalg.lstsq(A, dataY)[0]

T = np.matrix([T])

Ynew = V*T.T
Ynew = np.ravel(Ynew)


ax[2].plot(x, Ynew, 'r-')

plt.show()

