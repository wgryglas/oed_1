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

#import nlopt



#base = '/home/michal/avio/naca0012/single_sweapt/'
base = '/home/wgryglas/python/avio/naca0012/multi_sweapt_1/mach.20/'
#base = '/home/michal/avio/naca0012/multi_sweapt_1/all/'

fname0=base+'input/fin_%d.dat'
ofname0=base+'decomposed/fin_%d_modes_%d.dat'
oaggfname0=base+'aggregated/fin_%d_modes_%d.dat'
mfname0=base+'modes/mode_%d.dat'
mnpzfname0=base+'modes/mode_%d.npz'

bfile=base+'boundary.dat'
optimofile=base+'optim_data.npz'

experimentfile=base+'experiment.dat'
exp_aggfname=base+'reconstructed_experiment.dat'

NT = 10


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
experiment_data_raw = list()

for line in file(experimentfile):
   # line.split(" ")
    experiment_data_raw.append([float(xx) for xx in re.findall('[0-9,\.,\-,\+,e,E]+', line)])
        
               
 
experiment_data_raw = np.array(experiment_data_raw)[:,:]



def get_values_at_points(data,xx,yy):
    data_at_points = list()
    for _x, _y in zip(xx,yy):
        data_at_points.append( data[np.argmin((data[:,0] - _x)**2 + (data[:,1] - _y)**2), 2:] )

    return np.array(data_at_points)

#b_experiment = np.array(get_bs_for_xy(experiment_data_raw[:,0],experiment_data_raw[:,1]))


values_at_optim_points = get_values_at_points(experiment_data_raw,x_optim,y_optim)

#plt.plot(b_experiment, experiment_data_raw[:,3])
#plt.plot(x_optim, values_at_optim_points[:,1],'o-')
#plt.plot(experiment_data_raw[:,0], experiment_data_raw[:,3],'x-')


coeffs  = np.linalg.solve(A.T.dot(A), A.T.dot(np.array(values_at_optim_points[:,field_mesured])))
print coeffs
#coeffs2 = np.linalg.solve(A2.T.dot(A2), A2.T.dot( experiment_data[:,2+field_mesured] ))
for i in range(Nvars-1):
    f = np.zeros_like(all_vectors_b[0][:,field_mesured])
    
    plt.figure()
    plt.title('V%d'%i)
    plt.grid(which='both')
    for c,v in zip(coeffs, all_vectors_b):
        f = f + c *v[:,i]
      #  ff = ff + c2 *v[:,i]
   
    plt.plot(x_boundary, f,'o')
    plt.plot(experiment_data_raw[:,0], experiment_data_raw[:,2+i],'x')


############################################


modes = list()
modes_x = 0
modes_y = 0
for i in range(NT):
    mnpzfname = mnpzfname0%(i)
    
    fdata = np.load(mnpzfname)
    
    modes.append(fdata['v'])
    if i ==0:
        modes_x= fdata['x']
        modes_y= fdata['y']
doSave = True

spectrum_all = list()

deltaVi_all = list()
variables = ["V01", "V02", "V03", "V04", "V05", 'p']
#for j in range(1,NT+1):
#    print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
f =  np.zeros(N*Nvars)


for c,v in zip(coeffs, modes):
    v = np.array(v)

    f = f + c * v

if doSave:

    fp = file(exp_aggfname, 'w')
    fname = fname0%1       
    n = 0
    for line in file(fname):
        n = n + 1
        if n == 2:
            
            fp.write('VARIABLES = "X", "Y" ' + ",".join(['"'+_s+'_rec"' for _s in variables]) + '\n')
            
        elif n > 3 and  n < 4 + N:
            fp.write(" ".join([ str(modes_x[n-4]), str(modes_y[n-4]) ]))
            
            fp.write("    ")                
            
            #fp.write(" ".join([ str(x) for x in f.reshape(N,Nvars)[n-4,:]]))
            
            fp.write(" ".join([ str(q) for q in f.reshape(N,Nvars)[n-4,:] ]))
            
            fp.write("\n")
            
        else:
            fp.write(line)
            
    fp.close()
        



#plt.show()


print "DONE"