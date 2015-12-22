# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 15:43:16 2015

@author: mdzikowski
"""







import matplotlib as mpl

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
import modred as MR

from matplotlib.backends.backend_pdf import PdfPages


n = 0

AllData = list()

i = 1

#base = '/home/mdzikowski/avio/naca0012/single_sweapt'
#base = '/home/michal/avio/naca0012/single_sweapt/'
base = '/home/michal/avio/naca0012/multi_sweapt_1/mach.20/'
#base = '/home/michal/avio/naca0012/multi_sweapt_1/all/'
NT = 10

fname0=base+'input/fin_%d.dat'

ofname0=base+'decomposed/fin_%d_modes_%d.dat'
oaggfname0=base+'aggregated/fin_%d_modes_%d.dat'

mfname0=base+'modes/mode_%d.dat'

mnpzfname0=base+'modes/mode_%d.npz'

fname = fname0%i


pdf_pages = PdfPages(base+'/POD-stats.pdf')
pdf_nb_plots_per_page = 2
pdf_grid_size = (pdf_nb_plots_per_page, 1)


pltF = plt.figure
pltF = plt.figure
class plotToPdf:
    def __init__(self):
        self.i = 0

    def begin(self):
      if self.i % pdf_nb_plots_per_page == 0:
          self.fig = pltF(figsize=(8.27, 11.69), dpi=100)
     
      # Plot stuffs !
      plt.subplot2grid(pdf_grid_size, (self.i % pdf_nb_plots_per_page, 0))

      
    def end(self):
      if (self.i + 1) % pdf_nb_plots_per_page == 0:
        plt.tight_layout()
        pdf_pages.savefig(self.fig)        
      self.i = self.i + 1


pltM = plotToPdf()
pltM.begin()



def masker(**kargs):
    if len(kargs) == 0:
        pltM.end()
        pltM.begin()
    else:
        raise "Masker to much"        
plt.figure = masker

variables = list()
numpoints = 0
params = dict()
doSave = False

    
for line in file(fname):
    n = n + 1
    if re.match('.*VARIABLES.*', line):
        for var in re.findall('"([a-zA-Z,0-9]*)"',line):
            variables.append(var)
    else:
        for param in line.split(","):
            for m in re.findall('([a-zA-Z,0-9]*)=([a-zA-Z,0-9]*)',param)    :
                params[m[0]] = m[1]
    if n > 5:
        break
    n = 0 
N = int(params['N'])
R = 2
for i in range(1,NT+1):
    print i

    fname = fname0%i
    data = list()
    n = 0
    r = 0
    for line in file(fname):
        n = n + 1
        if n < 4:
            pass  
        else:
            #print line
            datarow = list()
            for num in re.findall('[0-9,\.,\-,e]+', line):
                datarow.append(float(num))
                
#            if (datarow[0]**2 + datarow[1]**2) < R**2:
            if True:
                data.append(datarow)
                
                    
            r = r + 1
        if r == N :
            print "done. ", len(data), " lines"
            break
        
    # add pressure
    kappa = 1.4        


    data = np.array(data)
    
    p = (kappa-1.) * data[:,2] * ( data[:,3] - (data[:,4]**2 + data[:,5]**2) / data[:,2] / 2. )
    data1 = np.zeros((N,len(data[0,:])+1))
    data1[:,:-1] = data
    data1[:,-1] = p
    
    AllData.append(data1)

variables.append('P')

N = len(AllData[0][:,0])


# 'flat' data









Nvars = len(AllData[0][0,:]) - 2

num_modes = NT


podInput = np.array([ rr[:,2:].reshape(N*Nvars) for rr in AllData ])
modes, eig_vals = MR.compute_POD_matrices_snaps_method(podInput.T, range(num_modes))


x = AllData[0][:,0]
y = AllData[0][:,1]




for i in range(num_modes):
    v = np.array(modes[:,i].T.tolist()[0])
    mnpzfname = mnpzfname0%(i)
    
    np.savez(mnpzfname, mode_num=i, N=N, Nvars=Nvars ,v=v, x=x, y=y)    
    
    fp = file(mfname0%(i), 'w')
    fname = fname0%int(1)        
    n = 0
    for line in file(fname):
        n = n + 1
        if n > 3 and  n < 4 + N:
            fp.write(" ".join([ str(q) for q in AllData[0][n-4,:2]]))
            fp.write("    ")                
            fp.write(" ".join([ str(q) for q in v.reshape(N,Nvars)[n-4,:]]))
            fp.write("\n")
        else:
            fp.write(line)




spectrum_all = list()

deltaVi_all = list()

pltM.begin()

for j in range(1,NT+1):
    print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    f =  np.zeros(N*Nvars)
    ff =  np.zeros(N*Nvars)    
    f0 = podInput[j-1,:]
    rhoE0 =  np.linalg.norm( f0.reshape(N,Nvars)[:,2] )
    energy = list()
    spectrum = list()
    
    deltaVi = [ list() for i in range(Nvars)]
    
    for i in range(num_modes):
        v = np.array(modes[:,i].T.tolist()[0])
        c = np.dot(v.T,f0)
        f = f + c * v
        ff = c * v
        #print j, i, np.linalg.norm(f-f0)
        
        rhoE = np.linalg.norm( f.reshape(N,Nvars)[:,2] )
        energy.append( rhoE / rhoE0 )
        spectrum.append( np.linalg.norm( (c*v).reshape(N,Nvars)[:,2] ) / rhoE0 )
        
        for l in range(Nvars):
            deltaVi[l].append( np.max( np.sqrt(np.power(f.reshape(N,Nvars)[:,l]- f0.reshape(N,Nvars)[:,l],2) ) ) )
                    
        
#  TITLE = ""
#VARIABLES = "X", "Y", "V01", "V02", "V03", "V04", "V05"
#ZONE T="", N=3704, E=7168, F=FEPOINT, ET=QUADRILATERAL
      
        if doSave:
            fp = file(ofname0%(j,i), 'w')
            ffp = file(oaggfname0%(j,i), 'w')
            fname = fname0%j        
            n = 0
            for line in file(fname):
                n = n + 1
                if n == 2:
                    
                    fp.write('VARIABLES = ' + ",".join(['"'+_s+'_over_'+_s+'_org"' for _s in variables]) + '\n')
                    ffp.write('VARIABLES = ' + ",".join(['"'+_s+'"' for _s in variables]) + '\n')
                elif n > 3 and  n < 4 + N:
                    fp.write(" ".join([ str(q) for q in AllData[0][n-4,:2]]))
                    ffp.write(" ".join([ str(q) for q in AllData[0][n-4,:2]]))
                    
                    fp.write("    ")                
                    ffp.write("    ")                
                    #fp.write(" ".join([ str(x) for x in f.reshape(N,Nvars)[n-4,:]]))
                    ffp.write(" ".join([ str(q) for q in ff.reshape(N,Nvars)[n-4,:] / f0.reshape(N,Nvars)[n-4,:] ]))
                    fp.write(" ".join([ str(q) for q in f.reshape(N,Nvars)[n-4,:] ]))
                    
                    fp.write("\n")
                    ffp.write("\n")
                else:
                    fp.write(line)
                    ffp.write(line)
            fp.close()
            ffp.close()
            
    deltaVi_all.append(deltaVi)                    
    spectrum_all.append(spectrum)
    plt.plot(energy,'o-')    

        #fname = fname0%[i,j]
plt.grid(which='both')    
plt.title(r'L2(V1)/L2(V1_orig)')
plt.xlabel('Number of modes used')
plt.ylabel('Percentage of energy restored')

plt.figure()
for s in spectrum_all:
    plt.plot(s,'o-')
plt.grid(which='both') 
plt.title(r'Energy carried by single mode')
plt.xlabel('Number of mode')
plt.ylabel('Relative ammount of energy')

for l in range(Nvars):    
    plt.figure()
    for m in range(NT):
        plt.plot(deltaVi_all[m][l],'o-')

    plt.title(r'Maximum diference for V%d'%l)
    plt.xlabel('Number of modes used')
    plt.ylabel('Maximum diference')
    plt.grid(which='both')    



#plt.show()        
pdf_pages.close()
    


print "DONE"