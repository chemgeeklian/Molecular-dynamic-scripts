# -*- coding: utf-8 -*-
"""
Created on Wed May  2 14:26:08 2018
ground and 1st excited state energy
@author: Xinran Lian
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
scaleu=35
scalev=26

files=list(np.empty(scaleu*scalev))
u=np.linspace(2.30,3.66,scaleu) #r1+r2
v=np.linspace(0.00,1.00,scalev) #r1-r2
stru=list(np.zeros(scaleu))
strv=list(np.zeros(scalev))
energy=np.zeros([scaleu,scalev,2]) #coordinate of u, v, and ground/ex1/ex2
transdip=np.zeros([scaleu,scalev,2]) #[:,:,0]: along x axis; [:,:,1]: y
u2d=np.empty([scaleu,scalev])
v2d=np.empty([scaleu,scalev])

for i in range(0,scaleu):
    stru[i]='%.2f'%u[i]
    for j in range(0,scalev):
        strv[j]='%.2f'%v[j]
        loc='groundlog/'
        files[i*scalev+j]=loc+"u"+stru[i]+"v"+strv[j]+"Au_grdp.log"

for i in range(0,scaleu):
    u2d[i,:]=u[i]
    for j in range(0,scalev):
        v2d[i,j]=v[j]
        rd=open(loc+"u"+stru[i]+"v"+strv[j]+"Au_grdp.log")
        for line in rd:
            if not line.find(' Ground to excited state transition electric dipole moments (Au):') == -1:
                line=rd.readline()
                line=rd.readline()
                transdip[i,j,0]=float(line[18:29])*2.541746
                transdip[i,j,1]=float(line[30:40])*2.541746
            if not line.find(' Excited State   1: ') == -1: #excited
                energy[i,j,1]=float(line[39:48])
            if not line.find('   N       SCF          CIS ') == -1: #ground
                line=rd.readline()
                line=rd.readline()
                energy[i,j,0]=float(line[8:20])
                continue
        rd.close()

energy[:,:,1]=energy[:,:,1]/27.2114+energy[:,:,0]

#make symmetry matrices, because r1 - r2 is symmetric left/right to 0
energysym= np.flip(energy[:,:,0],1)
energysym= np.append(energysym,energy[:,1:,0],axis=1).reshape([scaleu,-1])
energysym2=np.flip(energy[:,:,1],1) #excited energy
energysym2=np.append(energysym2,energy[:,1:,1],axis=1).reshape([scaleu,-1])
transdipx_sym=np.flip(transdip[:,:,0],1)
transdipx_sym=np.append(transdipx_sym,transdip[:,1:,0],axis=1).reshape([scaleu,-1])
transdipy_sym=np.flip(transdip[:,:,1],1)
transdipy_sym=np.append(transdipy_sym,transdip[:,1:,1],axis=1).reshape([scaleu,-1])

scalevsym=1+2*(scalev-1)

#make vector array for (r1-r2) in symmetric manner
vsym=np.flip(-v,0)
vsym=np.append(vsym,v[1:])
    
loc='surface/'
#Pot figure of energy surface
fig, ax = plt.subplots()
title='Transition dipole 01 sqrt(x^2+y^2)'
cax = ax.imshow(abs(transdipx_sym**2+transdipy_sym**2),extent=[-1*v[scalev-1], v[scalev-1], u[scaleu-1], 
                                  u[0]], interpolation='none', cmap=cm.coolwarm)
ax.set_xlabel('r1-r2 (Å)')
ax.set_ylabel('r1+r2 (Å)')
# Add colorbar, make sure to specify tick locations to match desired ticklabels
e01d=transdipx_sym.reshape(-1)
u1d=u2d.reshape(-1)
v1d=v2d.reshape(-1)
cbar = fig.colorbar(cax)
#plt.savefig(loc+title+'.png')
#ax.set_title(title+'\n')
plt.show()

#export matrix for further wavepacket calculation
'''np.savez('potentialenergy',pot1=energysym,pot2=energysym2)
np.savez('coordinate',u=u,v=vsym)
np.savez('transdipole',tsdpx=transdipx_sym,tsdpy=transdipy_sym)
'''