# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 17:35:21 2015

@author: michal
"""

class require:
    def __init__(self, *args, **kargs):
        for i,a in enumerate(args):
            if not a:
                raise "Rquired from-uppr argument number ", i , " failed"
        pass
    def __call__(self, func):
        return func
    

import matplotlib.pyplot as plt
import numpy as np
import bearded_octo_wookie.RED.RedReader as RedReader
import bearded_octo_wookie.RED.POD as POD
from scipy.spatial import cKDTree
from scipy.interpolate import interp1d 
from scipy.interpolate import UnivariateSpline


dir0 = '/home/michal/avio/naca0012/multi_sweapt_1/all/'
#dir0 = '/home/michal/avio/naca0012/multi_sweapt_1/mach.65/'
fname0 = dir0+'input/fin_%d.dat'
geometryFile = '/home/michal/avio/naca0012/multi_sweapt_1/all/name.get'
boundaryFile = '/home/michal/avio/naca0012/multi_sweapt_1/all/boundary.dat'
fpoints=dir0+'designed_points'
num_files = 120
num_modes = 15
kappa = 1.4


rtfs = [RedReader.RedTecplotFile(fname0%i) for i in range(1,num_files+1)]
    
#==============================================================================
# for rtf in rtfs:            
#     data = rtf.data        
#     p = (kappa-1.) * data[:,2] * ( data[:,3] - (data[:,4]**2 + data[:,5]**2) / data[:,2] / 2. )
#     rtf.appendData('p', p)          
#==============================================================================

modesH = POD.Modes(rtfs, num_modes=num_modes)    
modesH.writeModes('/tmp/test%d.vti')

mesh = cKDTree(modesH.baseRedFile.data[:,:2])

#==============================================================================
# SHOW RECONSTRUCTION 
#==============================================================================
plt.figure()
it = 0
modes_energy = np.zeros(num_modes)

numit = len(rtfs[0].variables[2:])*num_files*num_modes
for fid, field in enumerate(rtfs[0].variables[2:]):
    plt.subplot(3,2,fid+1)
    plt.title(field)
    en = np.zeros(num_modes)
    for i,rtf in enumerate(rtfs):
        shape = rtf.data[:,2:].shape
        shape1d = shape[0] * shape[1]
        v0 = rtf.data[:,2:].reshape((1,shape1d))
        v = np.zeros_like(v0)
        cs = list()
        print it ,'/' , numit
        
        for mid,m in enumerate(modesH.modes):
            
            it = it + 1
            c = v0.dot(m.reshape(shape1d).T)[0,0]
            v1 = c * m.reshape(shape1d)
            if fid ==0:
                en1 = np.sum( v1[:,1] + (v1[:,2]**2 + v1[:,3]**2 ) / v1[:,0] )
                en0 = np.sum( v0[:,1] + (v0[:,2]**2 + v0[:,3]**2 ) / v0[:,0] )
                modes_energy[mid] = modes_energy[mid] + en1/en0
            
            v = v + v1
        #    print v[:,fid]
            e0 = np.max(np.abs(v[:,fid] - v0[:,fid] ))
            en[mid] = en[mid] + e0
    plt.semilogy(en / num_files)
    plt.grid(which='both')            
        #print "EOF #############"
plt.figure()
plt.plot(modes_energy/ num_files, 'o')        
plt.grid(which='both')            
plt.show()