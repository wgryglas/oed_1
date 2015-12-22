# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 15:40:43 2015

@author: michal
"""

import numpy as N
import modred as MR
import matplotlib.pyplot as plt


num_vecs = 100

k = 0.1
# Arbitrary data
x,y = (1., 1.)
f = N.array((x,y))



vecs = [f + k * (0.5-N.random.random(2)) for i in range(num_vecs)]


podData = (N.array([ r.reshape(2) for r in vecs ])).T
num_modes = 2

print podData.shape

modes, eig_vals = MR.compute_POD_matrices_snaps_method(podData, range(num_modes))



f1 = N.zeros_like(f)



for m in modes:
    m = N.array(m)[0]
    c = m[0]*x + m[1]*y
    print c
    plt.plot([0.,c*m[0]],[0.,c*m[1]],'-')

plt.figure()

err = list()

for v in vecs:
    for m in modes:
        m = N.array(m)[0]
        v = N.array(v)
        c = m[0]*v[0] + m[1]*v[1]
        v = v - c * m
    err.append(v[0]*v[0]+v[1]*v[1])
plt.plot(err)
vecs = N.array(vecs)

plt.figure()
plt.plot(vecs[:,0],vecs[:,1],'o')
plt.plot(x,y,'o')

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