# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 15:40:43 2015

@author: michal
"""

import numpy as np
import modred as mr
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from wg.tools.function import *
from wg.tools.plot import *
#plt.show = empty

plt.figure = fixargument(plt.figure, figsize=(10, 10))

points = fixargument(plt.plot, override=False, color="red", linestyle="", marker="o", markersize=3)

def pntsDraw(x,y, *args):
    plt.figure()
    points(x,y, *args)
    latex_xticks("pi", np.pi, x, 0.01)
    plt.show()

pointsDraw = join(plt.figure, points, plt.show, argsFor=1)

num_vecs = 1000

k = 0.1
# Arbitrary data
x,y = (1., 1.)
f = np.array((x,y))


vecs = np.zeros((num_vecs, 2))

rand = fixargument(np.random.multivariate_normal, mean=[0, 0], cov=[[1, 0.99], [0.99, 1]])

#rand_norm = fixargument(np.random.multivariate_normal, mean=[0], cov=[[1]])
#rand = lambda : [rand_norm()[0], rand_norm()[0]]

vecs = np.array([rand() for i in range(num_vecs)])
#vecs[:, 1] = [rand()[1] for i in range(num_vecs)]



#vecs = [f + k * v for v in vecs]


#plt.figure()
#points([v[0] for v in vecs], [v[1] for v in vecs], color="red")


podData = (np.array([r.reshape(2) for r in vecs ])).T

print "pod data shape: ", podData.shape


num_modes = 2

modes, eig_vals = mr.compute_POD_matrices_snaps_method(podData, range(num_modes))

mods = np.array(modes)

print "modes shape", modes.shape
print "modes: ", modes
print "eig values shape", eig_vals.shape
print "eig values:", eig_vals



# res = np.zeros(2)
# for i in range(0, len(modes)):
#     res[i] = modes[i,0]*eig_vals[0] + modes[i,1]*eig_vals[1]
#
# print "results  ", res
#
# print "results 2", modes.dot(eig_vals)


plt.figure()
points(vecs[:,0], vecs[:,1], color="red")


@assert_arg_iterable(excepion=True, args=(0,1), kwargs=("X", "Y"))
def distance(X, Y):
    dist = np.zeros((len(X), len(X)))
    for i, (x1, y1) in enumerate(zip(X,Y)):
        for j,(x2, y2) in enumerate(zip(X,Y)):
            dist[i,j] = np.sqrt((x2-x1)**2 + (y2-y1)**2)

    if dist.shape[0] == 1 and dist.shape[1] == 1:
        return dist[0,0]
    else:
        return dist

#radius = max(distance(vecs[:, 0], vecs[:, 1]).reshape(vecs.shape[0]**2))/2
vecRange = join(fzip(max, min), substract)
rx = 0.9*vecRange(vecs[:, 0])/2
ry = 0.9*vecRange(vecs[:, 1])/2

cmap=plt.get_cmap('gist_rainbow')
for i in range(modes.shape[0]):
    m = modes[:, i]
    m = m.flat
    plt.arrow(0, 0, rx*m[0], ry*m[1], color=cmap(i*1./modes.shape[0]))

plt.show()

exit()

f1 = np.zeros_like(f)

#for i,m in enumerate(modes):





# for m in modes:
#     print "m=", m
#     m = np.array(m)[0]
#     print "m=", m, "\n--------"
#
#     c = m[0]*x + m[1]*y
#     print c
#     plt.plot([0., c*m[0]], [0., c*m[1]], '-')

plt.show()



plt.figure()

err = list()

for v in vecs:
    for m in modes:
        m = np.array(m)[0]
        v = np.array(v)
        c = m[0]*v[0] + m[1]*v[1]
        v = v - c * m
    err.append(v[0]*v[0]+v[1]*v[1])
plt.plot(err)
vecs = np.array(vecs)

plt.figure()
plt.plot(vecs[:,0], vecs[:,1], 'o')
plt.plot(x, y, 'o')

plt.show()



#A = modes
#print A.shape
#c = N.linalg.lstsq(A, f.reshape(nx*ny))[0]
#
#for i in range(num_modes):
#    v = modes[:,i]
#    f1 = f1 + c[i]*v.reshape(nx,ny)    
#    #print v.reshape(nx*ny).dot(f1.reshape(nx*ny,1))    
#    plt.figure()
#    plt.plot(f1.reshape(nx*ny).T,'-')
#    plt.plot(f.reshape(nx*ny))
#
#    plt.show()
#    #print N.max(f1)
#    #print N.min(f1)
#    #print "###################"
##    plt.figure()
##    res = N.power((f1-f).reshape(nx*ny),2)
##    print N.sqrt(N.sum(res))
##    plt.imshow(res.reshape(nx,ny))
##    plt.colorbar()
#plt.show()