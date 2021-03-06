# -*- coding: utf-8 -*-
"""
created on sat feb 28 15:40:43 2015

@author: michal
"""

import numpy as np
from matplotlib.pyplot import axes
import modred as mr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib2tikz import save as tikz_save

from wg.tools.function import *
from wg.tools.plot import *
#plt.show = empty

# plt.figure = fixargument(plt.figure, figsize=(5, 3.2))

points = fixargument(plt.plot, override=False, color="red", linestyle="", marker="o", markersize=5)

pointsdraw = join(plt.figure, points, plt.show, argsfor=1)


def compare_original_and_reduced(x1, x2, a):

    a.plot(x1[:,0], x1[:,1], color="red", markersize=2, marker="o", linestyle="")
    a.plot(x2[:,0], x2[:,1], color="blue", markersize=2, marker="o", linestyle="")
    # points(x1[:,0], x1[:,1], color="red", markersize=10)
    # points(x2[:,0], x2[:,1], color="blue", markersize=10)

    for x1, x2 in zip(x1, x2):
        dx = np.ravel(0.9*(x2-x1))
        x1 = np.ravel(x1)
        a.arrow(x1[0], x1[1], dx[0], dx[1], color='k', linewidth=1)
    a.grid(True)
    a.set_xlabel("$x_1$")
    a.set_ylabel("$x_2$")


# settings

num_vecs = 50
k = 0.1


# arbitrary data

sigx, sigy = (1, 1)
f = np.array((sigx, sigy))


rand = fixargument(np.random.multivariate_normal, mean=[0, 0], cov=[[sigx, 0.8], [0.8, sigy]])

#rand_norm = fixargument(np.random.multivariate_normal, mean=[0], cov=[[1]])
#rand = lambda : [rand_norm()[0], rand_norm()[0]]

vecs = np.array([rand() for i in range(num_vecs)])
#vecs[:, 1] = [rand()[1] for i in range(num_vecs)]



#vecs = [f + k * v for v in vecs]


#plt.figure()
#points([v[0] for v in vecs], [v[1] for v in vecs], color="red")


poddata = (np.array([r.reshape(2) for r in vecs ])).T

print "pod data shape: ", poddata.shape


num_modes = 2

modes, eig_vals = mr.compute_POD_matrices_snaps_method(poddata, range(num_modes))

modes = np.array(modes)

print "modes shape", modes.shape
print "modes: ", modes
print "eig values shape", eig_vals.shape
print "eig values:", eig_vals


d = np.matrix(np.zeros((2, 2)))
np.fill_diagonal(d, eig_vals)

dinv = np.matrix(np.zeros((2, 2)))
np.fill_diagonal(dinv, 1./eig_vals)

x = np.matrix(vecs.T)
u = x.T * modes * dinv
#
# print np.linalg.norm(modes*(d*u.t) - x)
#

b = np.array([modes[:, 0].dot(a) for a in poddata.T])


x2 = np.array([modes[:, 0]*c for c in b])

print x2

# x2 = (modes[:,0]*eig_vals[0]*u[:,0].T).T
#
#
fig = plt.figure()

axes2 = fig.add_subplot(1, 2, 1)
axes2.plot(vecs[:, 0], vecs[:, 1], color='b', linestyle="", marker=".", markersize=2)
xmean = np.mean(vecs[:, 0])
ymean = np.mean(vecs[:, 1])
axes2.arrow(xmean, ymean, modes[0, 0], modes[1, 0], color='k', linewidth=2, width=0.01)
axes2.arrow(xmean, ymean, modes[0, 1], modes[1, 1], color='k', linewidth=2, width=0.01)
axes2.grid()
axes2.set_xlabel("$x_1$")
axes2.set_ylabel("$x_2$")
axes2.set_aspect('equal')

axes1 = fig.add_subplot(1, 2, 2)
compare_original_and_reduced(vecs, x2, axes1)
axes1.set_aspect('equal')

fig.tight_layout()

tikz_save('/home/wgryglas/Desktop/prezentacja_esn/reduction_exmple.tex', figureheight='180', figurewidth='180')
plt.show()

exit()
#
#
# @assert_arg_iterable(excepion=true, args=(0,1), kwargs=("x", "y"))
# def distance(x, y):
#     dist = np.zeros((len(x), len(x)))
#     for i, (x1, y1) in enumerate(zip(x,y)):
#         for j,(x2, y2) in enumerate(zip(x,y)):
#             dist[i,j] = np.sqrt((x2-x1)**2 + (y2-y1)**2)
#
#     if dist.shape[0] == 1 and dist.shape[1] == 1:
#         return dist[0, 0]
#     else:
#         return dist
#
#radius = max(distance(vecs[:, 0], vecs[:, 1]).reshape(vecs.shape[0]**2))/2
vecrange = join(fzip(max, min), substract)
rx = 0.9*vecrange(vecs[:, 0])/2
ry = 0.9*vecrange(vecs[:, 1])/2

cmap=plt.get_cmap('gist_rainbow')
for i in range(modes.shape[0]):
    m = modes[:, i]
    plt.arrow(0, 0, rx*m[0], ry*m[1], color=cmap(i*1./modes.shape[0]))

plt.show()



# f1 = np.zeros_like(f)

# for m in modes:
#     print "m=", m
#     m = np.array(m)[0]
#     print "m=", m, "\n--------"
#
#     c = m[0]*x + m[1]*y
#     print c
#     plt.plot([0., c*m[0]], [0., c*m[1]], '-')

# plt.show()

# plt.figure()
#
# err = list()
#
# for v in vecs:
#     for m in modes:
#         m = np.array(m)[0]
#         v = np.array(v)
#         c = m.dot(v)#m[0]*v[0] + m[1]*v[1]
#         v = v - c * m
#     err.append(v[0]*v[0]+v[1]*v[1])
# plt.plot(err)
# vecs = np.array(vecs)
#
# plt.figure()
# plt.plot(vecs[:,0], vecs[:,1], 'o')
# plt.plot(x, y, 'o')
#
# plt.show()


x2 = np.array([rand() for i in range(num_modes)])

a = np.zeros((num_modes*2, num_modes*2))
a[0:num_modes, 0:num_modes] = modes
a[num_modes:2*num_modes, num_modes:2*num_modes] = modes

x2 = np.ravel(x2.reshape((1, 4)))

t = np.linalg.lstsq(a, x2)[0]
t = t.reshape((2, 2))

xnew = modes*t
xnew = xnew.T

print xnew
exit()


xc = np.matrix([[sigx], [sigy]])

c = np.linalg.lstsq(modes, xc)[0]

c = c.t


print "new coeffs", c


# reconstruction based on new data

dinv = np.matrix(np.zeros((num_modes, num_modes)))
np.fill_diagonal(dinv, 1./c)
d = np.matrix(np.zeros((num_modes, num_modes)))
np.fill_diagonal(d, c)

#u = x.t * modes * dinv
xnew = modes * d * u.t


print modes.t*xc*xc.t*modes

exit()

# print np.linalg.norm(x-xnew)

plt.figure()
points(x[0, :], x[1, :], color="blue", markersize=8)
points(xnew[0, :], xnew[1, :], color="red")
plt.show()

# a = modes
# print a.shape
# c = np.linalg.lstsq(a, f.reshape(nx*ny))[0]
#
# for i in range(num_modes):
#    v = modes[:,i]
#    f1 = f1 + c[i]*v.reshape(nx,ny)
#    #print v.reshape(nx*ny).dot(f1.reshape(nx*ny,1))
#    plt.figure()
#    plt.plot(f1.reshape(nx*ny).t,'-')
#    plt.plot(f.reshape(nx*ny))
#
#    plt.show()
#    #print n.max(f1)
#    #print n.min(f1)
#    #print "###################"
##    plt.figure()
##    res = n.power((f1-f).reshape(nx*ny),2)
##    print n.sqrt(n.sum(res))
##    plt.imshow(res.reshape(nx,ny))
##    plt.colorbar()
#plt.show()