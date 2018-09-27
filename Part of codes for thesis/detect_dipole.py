# -*- coding: utf-8 -*-
"""
Created on Thu May  3 11:50:06 2018
Permanent Dipole of ground and 1st excited state
@author: Xinran Lian
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
scaleu=35
scalev=26

u=np.linspace(2.30,3.66,scaleu) #r1+r2
v=np.linspace(0.00,1.00,scalev) #r1-r2
stru=list(np.zeros(scaleu))
strv=list(np.zeros(scalev))
dipole=np.zeros([scaleu,scalev,2]) #coordinate of u, v, and ground/ex1/ex2
u2d=np.empty([scaleu,scalev])
v2d=np.empty([scaleu,scalev])

locarray=['groundlog/','ex1log/']
for loc in (locarray):
    files=list(np.empty(scaleu*scalev))
    for i in range(0,scaleu):
        stru[i]='%.2f'%u[i]
        for j in range(0,scalev):
            strv[j]='%.2f'%v[j]
            files[i*scalev+j]=loc+"u"+stru[i]+"v"+strv[j]+"Au_grdp.log"
            
    for i in range(0,scaleu):
        u2d[i,:]=u[i]
        for j in range(0,scalev):
            v2d[i,j]=v[j]
            rd=open(loc+"u"+stru[i]+"v"+strv[j]+"Au_grdp.log")
            for line in rd:
                if not line.find(' Dipole moment (field-independent basis, Debye)') == -1: 
                    line=rd.readline()
                    dipole[i,j,locarray.index(loc)]=float(line[65:80])*2.541746
                    #convert dipole moment drom D to atomic unit
            rd.close()
            
dipolesym= np.flip(dipole[:,:,0],1)
dipolesym= np.append(dipolesym,dipole[:,1:,0],axis=1).reshape([scaleu,-1])
dipolesym2=np.flip(dipole[:,:,1],1) #excited energy
dipolesym2=np.append(dipolesym2,dipole[:,1:,1],axis=1).reshape([scaleu,-1])
scalevsym=1+2*(scalev-1)

#make vector array for (r1-r2) in symmetric manner
vsym=np.flip(-v,0)
vsym=np.append(vsym,v[1:])

loc='surface/'
title='Dipole moment (ground)'

fig, ax = plt.subplots(figsize=(10, 10*scaleu/scalevsym))
cax = ax.imshow(dipolesym,extent=[-1*v[scalev-1]*1.89, v[scalev-1]*1.89, u[scaleu-1]*1.89, 
                                  u[0]*1.89], interpolation='none', cmap=cm.coolwarm)
ax.set_title(title+'\n')
ax.set_xlabel('r1-r2 (a.u.)')
ax.set_ylabel('r1+r2 (a.u.)')
# Add colorbar, make sure to specify tick locations to match desired ticklabels
e01d=dipole[:,:,0].reshape(-1)
u1d=u2d.reshape(-1)
v1d=v2d.reshape(-1)
cbar = fig.colorbar(cax, ticks=[min(e01d), (min(e01d)+max(e01d))/2,max(e01d)])
cbar.ax.set_yticklabels([str(min(e01d)),
                         str((min(e01d)+round(max(e01d),5))/2), round(max(e01d),5)])
plt.savefig(loc+title+'.png')
plt.show()

#export matrix for further wavepacket calculation
#np.savez('dipole',dip1=dipolesym,dip2=dipolesym2)