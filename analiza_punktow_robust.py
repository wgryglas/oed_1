# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 15:37:22 2015

@author: mdzikowski
"""

import numpy as np
import matplotlib.pyplot as plt

import bearded_octo_wookie.RED.RedReader as RedReader


from scipy.spatial import cKDTree
from scipy.interpolate import interp1d 
from scipy.interpolate import UnivariateSpline


geometryFile = '/home/michal/avio/naca0012/multi_sweapt_1/all/name.get'
boundaryFile = '/home/michal/avio/naca0012/multi_sweapt_1/all/boundary.dat'
#exerimentFile = '/home/mdzikowski/avio/naca0012/multi_sweapt_1/aoa2.2_ma_0_2.dat'
#exerimentFile = '/home/mdzikowski/avio/naca0012/multi_sweapt_1/aoa2_2_ma_0_7.dat'
exerimentFile = '/home/mdzikowski/avio/naca0012/multi_sweapt_1/mach.65/input/fin_1.dat'
#fpoints = '/home/mdzikowski/avio/naca0012/multi_sweapt_1/mach.20/designed_points.npz'
fpoints = '/home/mdzikowski/avio/naca0012/multi_sweapt_1/all/designed_points.npz'
#fpoints = '/home/mdzikowski/avio/naca0012/multi_sweapt_1/mach.65/designed_points.npz'
#fpoints = '/home/mdzikowski/avio/naca0012/multi_sweapt_1/all/designed_points.npz'


work_on_param = 1

get = RedReader.GetReader(geometryFile)

boundary_nodes = np.loadtxt(boundaryFile)[:,:2]
boundary_nodes_params = np.array(get.getParamList(boundary_nodes[:,:2].tolist(), 1))


experiment_all = RedReader.RedTecplotFile(exerimentFile) 

mesh = cKDTree(experiment_all.data[:,:2])

bids = mesh.query(boundary_nodes)[1]    


boundary_data = experiment_all.data[bids,:]

boundary_data_param = get.getParamList(boundary_data[:,:2].tolist(),1)
boundary_values = boundary_data[:,2+work_on_param]

experiment_interpolator = interp1d(boundary_data_param, boundary_values, kind='cubic')

measure_points_file = np.load(fpoints)
measures_points = measure_points_file['testpoints_all_combinations']
information_matrixes = measure_points_file['information_matrixes']
modes = np.array(measure_points_file['boundary_modes'])
modes_points = np.array(measure_points_file['modes_points'])

measure_points = measures_points[-1]
A = information_matrixes[-1]

measure_points_uniform = np.linspace(0.,1., len(measures_points[0]))
A_uniform = information_matrixes[0]



realizations = list()
realizations_uniform = list()
realizations_polinom = list()
realizations_base = list()
realizations_polinom_uniform = list()


epsilon0 = 0.05
#epsilon0 = epsilon0*(np.max(boundary_values) - np.min(boundary_values))

plt.figure()
for m in modes:
    plt.plot(modes_points, m)
plt.figure()


for i in range(100):
    #delta = epsilon0 * (np.random.random(len(boundary_values)) - 0.5)
    delta = epsilon0 * (np.random.normal(size=len(measure_points))  )
    delta[:3]=0
    delta[-3:]=0
    experiment_interpolator = interp1d(boundary_data_param, boundary_values, kind='cubic',bounds_error=True )
    

    tested_values = experiment_interpolator(measure_points+delta)
    tested_values_uniform = experiment_interpolator(measure_points_uniform+delta)
    
    cs = np.linalg.solve(A.T.dot(A),A.T.dot(tested_values))
    realizations.append(cs.dot(modes))

    cs_uniform = np.linalg.solve(A_uniform.T.dot(A_uniform),A_uniform.T.dot(tested_values_uniform))
    realizations_uniform.append(cs_uniform.dot(modes))

    #realizations_polinom.append(interp1d(measure_points, tested_values, bounds_error=False, kind='scubi')(modes_points))

    realizations_polinom_uniform.append(interp1d(measure_points_uniform, tested_values_uniform, bounds_error=False, kind='linear')(modes_points))
    
    
    realizations_base.append(experiment_interpolator(modes_points))
    #realizations_polinom.append(experiment_interpolator(modes_points))
 #   plt.plot(measure_points, tested_values, 'o')
    #plt.plot(modes_points, cs.dot(modes), '-', lw=1)
fig, (ax) = plt.subplots(1, 1, sharex=True)

#==============================================================================

x = modes_points
y0 = np.mean(realizations_base, axis=0)
stddervi = np.std(realizations_base, axis=0, ddof=1)
y1 = y0 - stddervi
y2 = y0 + stddervi

ax.plot(boundary_nodes_params,boundary_values,'b-', label='"Real"')
ax.plot(x,y1 , 'b--',lw=1)
ax.plot(x,y2 , 'b--',lw=1)
ax.fill_between(x, y1, y2, where=y2>=y1, facecolor='b', interpolate=True, alpha=0.2)
#==============================================================================

x = modes_points
y0 = np.mean(realizations, axis=0)
stddervi = np.std(realizations, axis=0, ddof=1)
y1 = y0 - stddervi
y2 = y0 + stddervi

ax.plot(x,y0 ,'ko', label='Reconstructed - optimal')
ax.plot(x,y1 , 'k--',lw=1)
ax.plot(x,y2 , 'k--',lw=1)
ax.fill_between(x, y1, y2, where=y2>=y1, facecolor='k', interpolate=True, alpha=0.2)

#==============================================================================

x = modes_points
y0 = np.mean(realizations_uniform, axis=0)
stddervi = np.std(realizations_uniform, axis=0, ddof=1)
y1 = y0 - stddervi
y2 = y0 + stddervi

ax.plot(x,y0 ,'wo', label='Reconstructed - uniform')
ax.plot(x,y1 , 'r--',lw=1)
ax.plot(x,y2 , 'r--',lw=1)
ax.fill_between(x, y1, y2, where=y2>=y1, facecolor='r', interpolate=True, alpha=0.2)

#==============================================================================

#==============================================================================

x = modes_points
y0 = np.mean(realizations_polinom_uniform, axis=0)
stddervi = np.std(realizations_polinom_uniform, axis=0, ddof=1)
y1 = y0 - stddervi
y2 = y0 + stddervi

ax.plot(x,y0 ,'g-', label='Polinominal - uniform')
ax.plot(x,y1 , 'g--',lw=1)
ax.plot(x,y2 , 'g--',lw=1)
ax.fill_between(x, y1, y2, where=y2>=y1, facecolor='g', interpolate=True, alpha=0.2)


zero = np.average(realizations)
for p in measure_points:
    plt.plot((p,p), (-1000, zero), 'k--', lw=1)
for p in measure_points_uniform:
    plt.plot((p,p), (zero, 1000), 'k--', lw=1)
#==============================================================================
# plt.plot(modes_points, np.mean(realizations_polinom, axis=0),'c', label='Interpolated - optimal')
# plt.plot(modes_points, np.min(realizations_polinom, axis=0), 'c--',lw=1)
# plt.plot(modes_points, np.max(realizations_polinom, axis=0), 'c--',lw=1)
#==============================================================================




plt.ylim(np.min(boundary_values - 0.7 * epsilon0), np.max(boundary_values + 0.7 * epsilon0))
plt.legend()
plt.xlabel('Position on surface, 0.5=L.E')
plt.ylabel('Nondimensional total pressure')

plt.figure()


rec_eps = np.std(realizations, axis=0, ddof=1)
plt.plot(modes_points,rec_eps, label='Reconstructed - optimal')


rec_eps = np.std(realizations_uniform, axis=0, ddof=1)
plt.plot(modes_points,rec_eps, label='Reconstructed - uniform')

rec_eps = np.std(realizations_polinom_uniform, axis=0, ddof=1)
plt.plot(modes_points,rec_eps, label='Polinominal - uniform')

plt.ylabel('Std. dervi.')
plt.xlabel('Position on surface, 0.5=L.E')

rec_eps_base = np.std(realizations_base, axis=0, ddof=1)
plt.plot(modes_points, rec_eps_base, label='Base')
plt.legend()

#==============================================================================
# rec_avg_uniform = np.mean(realizations, axis=0)
# rec_eps_uniform = np.average(realizations_uniform - rec_avg_uniform, axis=0)
# plt.plot(rec_eps_uniform, label='Reconstructed - uniform')
#==============================================================================


#==============================================================================
# rec_avg_polinom_uniform = np.average(realizations_polinom_uniform, axis=0)
# rec_min_polinom_uniform = np.min(realizations_polinom_uniform, axis=0)
# rec_max_polinom_uniform = np.max(realizations_polinom_uniform, axis=0)
# rec_eps_polinom_uniform = (rec_max_polinom_uniform - rec_min_polinom_uniform) 
# plt.plot(rec_eps_polinom_uniform, label='Interpolated - uniform')
#==============================================================================

#==============================================================================
# rec_avg_polinom = np.average(realizations_polinom, axis=0)
# rec_min_polinom = np.min(realizations_polinom, axis=0)
# rec_max_polinom = np.max(realizations_polinom, axis=0)
# rec_eps_polinom = (rec_max_polinom - rec_min_polinom) / rec_avg_polinom
# plt.plot(rec_eps_polinom, label='Interpolated - optimal')
#==============================================================================



#plt.ylim( -0.5 * epsilon0,  0.5 * epsilon0)





plt.show()
