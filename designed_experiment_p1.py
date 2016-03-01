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



#base = '/home/wgryglas/python/avio/naca0012/single_sweapt'
#base = '/home/wgryglas/python/avio/naca0012/single_sweapt/'
base = '/home/wgryglas/python/avio/naca0012/multi_sweapt_1/mach.20/'
NT = 10

#base = '/home/wgryglas/python/avio/naca0012/multi_sweapt_1/all/'
#NT = 120

num_modes = 3

field_mesured = 1
num_test_points = 64


fname0=base+'input/fin_%d.dat'
ofname0=base+'decomposed/fin_%d_modes_%d.dat'
oaggfname0=base+'aggregated/fin_%d_modes_%d.dat'
mfname0=base+'modes/mode_%d.dat'
mnpzfname0=base+'modes/mode_%d.npz'

bfile=base+'boundary.dat'
optimofile=base+'optim_data.npz'


vectors = list()

all_vectors = list()

mnpzfname = mnpzfname0%(0)
fdata = np.load(mnpzfname)

N = fdata['N']
Nvars = fdata['Nvars']
x = np.array(fdata['x'])
y = np.array(fdata['y'])



for i in range(NT):
    mnpzfname = mnpzfname0%(i)
    
    fdata = np.load(mnpzfname)
    
    all_vectors.append(fdata['v'])
    
 
 
#read boundary nodes list
bnodes = list()

for line in file(bfile):
   # line.split(" ")
    xn, yn, idn = re.findall('[0-9,\.,\-,\+,e,E]+', line)
    idn = int(float(idn)) - 1
    
    bnodes.append(idn)

    #check!!!!
    if (float(xn)-x[idn])**2 + (float(yn)-y[idn])**2 > 1E-10:
        print "ERROR, node: ", idn 
               
 
 
bnodes = np.array(bnodes)



print "DONE"
xn = x[bnodes]
yn = y[bnodes]

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
bnodes1 = bnodes[sorted_keys]

shit_to_cont = -np.argmax(rn)

xn1 = np.roll(xn1, shit_to_cont)
yn1 = np.roll(yn1, shit_to_cont)
rn = np.roll(rn, shit_to_cont)
alpha = np.roll(alpha, shit_to_cont)
bnodes1 = np.roll(bnodes1,shit_to_cont)

l = 0.

bparam = list()
bparam.append(0.)
for i in range(1, len(bnodes)):
    dl = np.sqrt( (xn1[i-1]-xn1[i])**2 + (yn1[i-1]-yn1[i])**2  )
    l = l +dl
    bparam.append(l)


dl = np.sqrt( (xn1[-1]-xn1[0])**2 + (yn1[-1]-yn1[0])**2  )
l = l +dl
bparam.append(l)



bparam = np.array(bparam)
bparam = bparam / np.max(bparam)
bparam[-1] = 1.
bparam1 = bparam[:-1]
#spline representation

tck_x = interpolate.splrep(bparam, np.hstack([xn1,xn1[0]]), s=0.)
tck_y = interpolate.splrep(bparam, np.hstack([yn1,yn1[0]]), s=0.)

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
    
#plt.figure()
#
#
#plt.plot(x_new, y_new, '-')
#
#plt.plot(xn,yn, 'x')
#
#plt.show()        



all_vectors_b = list()
all_vectors_b_interp = list()
for v in all_vectors:
    v1 = v.reshape(N,Nvars)[bnodes1,:]
    all_vectors_b.append(v1)

#    fs = list()
#    for j in range(Nvars):
#        fs.append( interpolate.interp1d(bparam1, v1[:,j]) )
#    all_vectors_b_interp.append(fs)
    
    all_vectors_b_interp.append( interpolate.interp1d(bparam1.T, v1.T) )
    
#for i in range(num_modes):
#    plt.plot(bparam1, all_vectors_b[i][:,0],'o')        
#    plt.plot(bparam_new, all_vectors_b_interp[i](bparam_new)[0,:])        

#for i in range(Nvars):
#    f = np.zeros_like(all_vectors_b[i][:,0])
#    for m,v in enumerate(all_vectors_b):
#        c = v[:,i].dot(all_vectors_b[m][:,i])
#        f = f + c * v[:,i]
#    plt.figure()
#    plt.plot(bparam1,f,'o')
#    plt.plot()
#    plt.plot(bparam1, all_vectors_b[i][:,0],'o')   
#    plt.plot(bparam1, all_vectors_b[i][:,0],'o')        
#   plt.plot(bparam_new, all_vectors_b_interp[i](bparam_new)[0,:])        

#plt.show()
##picking
    
    




def opt_to_bparams(_opt_params):
    _bparams = list()
    p0 = 0.
    for p in _opt_params:
        p0=p/2.+p0
        _bparams.append(p0)
        p0=p/2.+p0
    return 0.99 * np.array(_bparams) / p0 
    
def information_matrix(_bparams):
    _A = np.zeros((num_test_points,num_modes))

    for i in range(num_modes):
        for j in range(num_test_points):
            _A[j,i] = all_vectors_b_interp[i](_bparams[j])[field_mesured]       
    return _A
def target_f(opt_params,grad):

    _A = information_matrix(opt_to_bparams(opt_params))
    _A = _A.T.dot(_A)
#    if not num_test_points == num_modes:
#        _A =  _A.T.dot(_A)
        
    try:
        r = np.trace(np.linalg.inv(_A))
        print r
        return r
    except np.linalg.linalg.LinAlgError:
        return np.nan
    

def target_f_nograd(bparams):
    _A = information_matrix(bparams)
    try:
        r = np.trace(np.linalg.inv(_A))
        print r*r
        return r*r
    except np.linalg.linalg.LinAlgError:
        return np.nan

    
#bs = 1.-np.linspace(0.0011,0.99,num_test_points)
bs = np.ones(num_test_points)

print opt_to_bparams(bs)


   


opt = nlopt.opt(nlopt.LN_COBYLA, num_test_points)
#opt = nlopt.opt(nlopt.LD_CCSAQ, num_test_points)

opt.set_lower_bounds(np.ones(num_test_points)*0.001)
opt.set_upper_bounds(2.*np.ones(num_test_points))
opt.set_min_objective(target_f)

#opt.inequality_constraint(myconstraint1, 1e-6)
#opt.add_inequality_constraint(myconstraint2, 1e-6)
#opt.add_equality_constraint(myconstraint1,1E-4)

opt.set_xtol_rel(1e-3)
opt.set_maxeval(1)

bs_optim = opt_to_bparams( opt.optimize(bs.tolist()) )
#bs_optim = bparam1
minf = opt.last_optimum_value()

print "minimum value = ", minf
print "result code = ", opt.last_optimize_result()
    
plt.figure()


plt.plot(x_new, y_new, '-')

x_optim = interpolate.splev(bs_optim, tck_x, der=0) + x0
y_optim = interpolate.splev(bs_optim, tck_y, der=0) + y0

plt.plot(x_optim, y_optim, 'o', label='Optimal points')

#plt.plot(xn,yn, 'x')



x_bs = interpolate.splev(opt_to_bparams(bs), tck_x, der=0) + x0
y_bs = interpolate.splev(opt_to_bparams(bs), tck_y, der=0) + y0

plt.plot(x_bs, y_bs, '^',label='Starting points')
plt.legend()
plt.grid(which='both')


np.savez(optimofile, 
#        A=information_matrix(opt_to_bparams(bs_optim)), 
        A=information_matrix(bs_optim), 
        x_optim=x_optim, 
        y_optim=y_optim, 
        field_mesured=field_mesured, 
        num_modes=num_modes,
        num_test_points=num_test_points,
        A_opt_start=information_matrix(opt_to_bparams(bs)),
        all_vectors_b=all_vectors_b,
        x_boundary = xn1+x0,
        y_boundary = yn1+y0
)

x_boundary = xn1+x0
y_boundary = yn1+y0
A=information_matrix(bs_optim)

plt.figure()





###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################




#Aorg = information_matrix(bs)
#
#A = A.T.dot(A)
#Aorg = Aorg.T.dot(Aorg)
#plt.figure()
#
#eAorg = np.linalg.eigvals(Aorg)
#eA = np.linalg.eigvals(A)
#
#
#
#plt.semilogy(np.real(eA),'o', label='Eigenvals optimal')
#plt.semilogy(np.real(eAorg),'x', label='Eigenvals start')
#
#plt.legend()
#
#plt.grid(which='both')
#
#
#plt.figure()
#
#plt.plot(bs_optim, 'o', label='Optimal points')
#plt.plot(bs, '^',label='Starting points')
#plt.grid(which='both')
#plt.legend()
#
#plt.show()
#
















#    


#A =  target_f(bs)
#
#cs = np.linspace(-1.,1., num_modes)
#
#plt.plot(bs, A.dot(cs),'o')
#
#
#f = np.zeros_like(all_vectors[0])
#
#for c,v in zip(cs,all_vectors):
#    f = f + c * v
#
#plt.plot(bparam1, f.reshape(N,Nvars)[bnodes1,mode_mesured],'-')
#

#
#    
#    
plt.show()

print "DONE"