# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 15:40:43 2015

@author: michal
"""

import numpy as N
import modred as MR
import matplotlib.pyplot as plt

nx = 200
ny = 200
num_vecs = 10
dlen = nx*ny
k = 0.1
# Arbitrary data
x,y = N.meshgrid(N.linspace(0,1.,ny), N.linspace(0,1.,nx))



f = N.sin(2. * N.pi * y)  * N.sin(2. * N.pi * x) 
vecs = [f + k * N.random.random((nx, ny)) for i in range(num_vecs)]


podData = (N.array([ r.reshape(nx*ny) for r in vecs ])).T
num_modes = num_vecs
#print vecs.shape
print podData.shape
modes, eig_vals = MR.compute_POD_matrices_snaps_method(podData, range(num_modes))


#plt.imshow(f)
#plt.figure()
#plt.imshow(vecs[0])
#
#plt.figure()
#plt.plot(podData[0])
#plt.plot(f.reshape(nx*ny))
#plt.show()
f1 = N.zeros_like(f)
print eig_vals
print modes.shape


A = modes
print A.shape
c = N.linalg.lstsq(A, f.reshape(nx*ny))[0]

for i in range(num_modes):
    v = modes[:,i]
    f1 = f1 + c[i]*v.reshape(nx,ny)    
    #print v.reshape(nx*ny).dot(f1.reshape(nx*ny,1))    
    plt.figure()
    plt.plot(f1.reshape(nx*ny).T,'-')
    plt.plot(f.reshape(nx*ny))

    plt.show()
    #print N.max(f1)
    #print N.min(f1)
    #print "###################"
#    plt.figure()
#    res = N.power((f1-f).reshape(nx*ny),2)
#    print N.sqrt(N.sum(res))
#    plt.imshow(res.reshape(nx,ny))
#    plt.colorbar()
plt.show()