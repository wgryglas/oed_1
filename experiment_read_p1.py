# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 15:43:16 2015

@author: mdzikowski
"""







import matplotlib as mpl
from scipy import interpolate
from scipy import optimize

mpl.rcParams['lines.linewidth'] = 4
mpl.rcParams['lines.markersize'] = 12
mpl.rcParams['lines.markeredgewidth'] = 2.
mpl.rcParams['lines.color'] = 'w'
font = {
        'size'   : 20}

mpl.rc('font', **font)

import numpy as np
import matplotlib.pyplot as plt
import re

import nlopt



base = '/home/michal/avio/naca0012/single_sweapt/'

fname0=base+'input/fin_%d.dat'
ofname0=base+'decomposed/fin_%d_modes_%d.dat'
oaggfname0=base+'aggregated/fin_%d_modes_%d.dat'
mfname0=base+'modes/mode_%d.dat'
mnpzfname0=base+'modes/mode_%d.npz'

bfile=base+'boundary.dat'
optimofile=base+'optim_data.npz'

experimentfile=base+'experiment.dat'


NT = 20


vectors = list()

all_vectors = list()

mnpzfname = mnpzfname0%(0)
fdata = np.load(mnpzfname)

N = fdata['N']
Nvars = fdata['Nvars']

optdata = np.load(optimofile)

A=optdata['A']
x_optim=optdata['x_optim']
y_optim=optdata['y_optim']
field_mesured=optdata['field_mesured']
num_modes=optdata['num_modes']
num_test_points=optdata['num_test_points']
all_vectors_b=optdata['all_vectors_b']
x_boundary = optdata['x_boundary']
y_boundary = optdata['y_boundary']


for i in range(NT):
    mnpzfname = mnpzfname0%(i)
    
    fdata = np.load(mnpzfname)
    
    all_vectors.append(fdata['v'])
    
 
 
#read boundary nodes list and data
experiment_data = list()

for line in file(experimentfile):
   # line.split(" ")
    experiment_data.append([float(xx) for xx in re.findall('[0-9,\.,\-,\+,e,E]+', line)])
        
               
 
experiment_data = np.array(experiment_data)



xn = experiment_data[:,0]
yn = experiment_data[:,1]
bnodes = np.arange(0,len(xn))


x0 = np.sum(xn) / len(xn)
y0 = np.sum(yn) / len(xn)

xn1 = xn - x0
yn1 = yn - y0

rn = np.sqrt(yn1**2 + xn1**2)
alpha = 1./np.arctan2(xn1,yn1)


sorted_keys = np.argsort(alpha)

xn1 = xn1[sorted_keys]
yn1 = yn1[sorted_keys]
rn = rn[sorted_keys]
alpha = alpha[sorted_keys]
bnodes = bnodes[sorted_keys]



shit_to_cont = -np.argmax(rn)

xn1 = np.roll(xn1, shit_to_cont)
yn1 = np.roll(yn1, shit_to_cont)
rn = np.roll(rn, shit_to_cont)
alpha = np.roll(alpha, shit_to_cont)
bnodes = np.roll(bnodes, shit_to_cont)




l = 0.

bparam = list()
bparam.append(0.)



for i in range(1,len(xn1)):
    dl = np.sqrt( (xn1[i-1]-xn1[i])**2 + (yn1[i-1]-yn1[i])**2  )
    l = l +dl
    bparam.append(l)


dl = np.sqrt( (xn1[-1]-xn1[0])**2 + (yn1[-1]-yn1[0])**2  )
l = l + dl
bparam.append(l)



bparam = np.array(bparam)
bparam = bparam / np.max(bparam)
bparam[-1] = 1.
bparam1 = bparam[:-1]
#spline representation

tck_x = interpolate.splrep(bparam, np.hstack([xn1,xn1[0]]), s=0.)
tck_y = interpolate.splrep(bparam, np.hstack([yn1,yn1[0]]), s=0.)




boundary_interpolate = interpolate.interp1d(bparam1.T, experiment_data[bnodes,2:].T,kind='cubic') 

bparam_new = np.linspace(0.01,0.99,1000.)
x_new = interpolate.splev(bparam_new, tck_x, der=0) + x0
y_new = interpolate.splev(bparam_new, tck_y, der=0) + y0    



def get_bs_for_xy(xx,yy):

    #probe space, to have good starting point
    bs = np.linspace(0.,1.,num_test_points)
    br = list()
    for xs,ys in zip(xx,yy):
        f = lambda(b):  (interpolate.splev(b, tck_x, der=0) + x0  - xs)**2 + (interpolate.splev(b, tck_y, der=0) + y0  -ys)**2   
        #b_optim.append(optimize.fmin_cobyla(f, [bs[np.argmin(f(bs))]], [lambda(_x): _x>0., lambda(_x): _x< 1.], rhoend=1e-16)[0])
        br.append(optimize.fmin(f, bs[np.argmin(f(bs))])[0])
    return br        
        
b_optim = np.array(get_bs_for_xy(x_optim, y_optim)).T
b_boundary = np.array(get_bs_for_xy(x_boundary, y_boundary)).T


#plt.plot(bparam1,  experiment_data[bnodes,1], 'o')
#plt.plot(bparam_new,  boundary_interpolate(bparam_new)[1,:], 'x')
#plt.plot(b_optim, y_optim)
#plt.show()

#plt.plot(bparam)
#plt.plot(get_bs_for_xy(xn1+x0,yn1+y0),'o')
#plt.show()

#plt.plot(b_optim)
#plt.plot(b_boundary,'o')
#plt.plot(bparam1)
#plt.show()

#plt.figure()

print b_optim
val_at_test_points = boundary_interpolate(b_optim)[field_mesured,:]
#val_at_test_points =  experiment_data[bnodes,2+field_mesured]

#plt.plot(x_optim, 'o')
#plt.plot(bparam_new, boundary_interpolate(bparam_new)[field_mesured,:], '-')
#plt.plot(x_boundary, 'x')
#plt.show()


#print A.T.dot(A), A.T.dot(val_at_test_points)




A2 = list()

for m,v in enumerate(all_vectors_b):
    A2.append( v[:,field_mesured] )


A2 = np.array(A2).T

plt.figure()

coeffs  = np.linalg.solve(A.T.dot(A),    A.T.dot(val_at_test_points)                  )
#coeffs2 = np.linalg.solve(A2.T.dot(A2), A2.T.dot( experiment_data[:,2+field_mesured] ))


for i in range(Nvars-1):
    f = np.zeros_like(all_vectors_b[0][:,field_mesured])
    ff = np.zeros_like(all_vectors_b[0][:,field_mesured])
    plt.figure()
    plt.title('V%d'%i)
    plt.grid(which='both')
    for c,v in zip(coeffs, all_vectors_b):
        f = f + c *v[:,i]
      #  ff = ff + c2 *v[:,i]
   
    plt.plot(bparam1, experiment_data[:,2+i], 'o')
    plt.plot(b_boundary, f, '^')
    #plt.plot(bparam1, ff, 'x')
#plt.figure()



#plt.plot(x_boundary,y_boundary,'o')
#plt.plot(experiment_data[:,0], experiment_data[:,1],'x')





f = np.zeros_like(all_vectors_b[0])

for v , coeff in zip(coeffs, all_vectors_b):
    f = f + v * coeff

#f = f.reshape(Nvars,len(all_vectors_b[0]) / Nvars)

#for i in range(Nvars-1):
#    plt.figure()
#
#    plt.plot(bparam_new, boundary_interpolate(bparam_new)[i,:], '-')
#    plt.plot(b_optim, boundary_interpolate(b_optim)[i,:], 'o')
#    plt.plot( get_bs_for_xy(x_boundary, y_boundary), f[:,i],'*' )
    
plt.show()


print "DONE"