# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 20:22:39 2018

@author: Xinran Lian
"""

import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
from matplotlib import cm

zip=np.load('coordinate.npz')
energy_load=np.load('potentialenergy.npz')
energy=energy_load['pot2']
u=np.round(zip['u']*1.8897261339213,2) #r1 + r2
v=np.round(zip['v']*1.8897261339213,2) #r1 - r2 convert armstrong to atomic unit
du=round(u[2]-u[1],2)
dv=round(v[2]-v[1],2)
miu1= 1836 #4m(H)
miu2= 1836*(9*1836)/(11*1836)
scaleu=np.size(energy[:,0])
scalev=np.size(energy[0])

eyeu=sparse.eye(scalev)
eyev=sparse.eye(scaleu)

pot=sparse.eye(scaleu*scalev).toarray()*(energy.T).reshape(-1)

kiu=sparse.diags([1, -2, 1], [-1, 0, 1], shape=(scaleu, scaleu)).todense()
kiv=sparse.diags([1, -2, 1], [-1, 0, 1], shape=(scalev, scalev)).todense()
kiu_kron=sparse.kron(eyeu,kiu).toarray() #kinetic matrix
kiv_kron=sparse.kron(kiv,eyev).toarray()

totalkiu=(-1/miu1)*kiu_kron/(du**2)
totalkiv=(-1/miu2)*kiv_kron/(dv**2)
ham=totalkiu+totalkiv+pot  #Hamitonian
eign,vector=np.linalg.eig(ham) 
'''such that the column vector[:,i] is the eigenvector corresponding to the 
eigenvalue eign[i]

eignandvector:sort vectors according to value of eignvalue from low to high
'''
eignandvector=np.append(eign,vector).reshape([-1,scaleu*scalev])
eignandvector_sort=eignandvector.T[eignandvector.T[:,0].argsort(kind='mergesort')].T
number=3
e1=eignandvector_sort[1:,number].reshape([-1,35])
e1=e1.T 

#check number of nodes (ground no, 1st excited 1 node...)
#vibrational excited states.
#reshape each column ov vector

#print(np.all(ham.T==ham))
eign=np.sort(eign,kind='mergesort')
if np.all(energy==energy_load['pot1']):
    state_str='ground'
else:
    if np.all(energy==energy_load['pot2']):
        state_str='1st excited'
    else:
        state_str='[ERROR]'

np.save(state_str+'_eigens',eignandvector_sort[:,:20]) #enough to save first 20 eigen states

e=e1.reshape(-1)
fig, ax = plt.subplots(figsize=(8, 8*35/51))
cax = ax.imshow(e1,extent=[-1*v[scalev-1], v[scalev-1], u[scaleu-1], u[0]], 
                interpolation='nearest', cmap=cm.coolwarm)
ax.set_title('Wavefunction of '+state_str+' BeH2 No. '+str(number)+' vector\n')
ax.set_xlabel('r1-r2 (a.u.)')
ax.set_ylabel('r1+r2 (a.u.)')
# Add colorbar, make sure to specify tick locations to match desired ticklabels
cbar = fig.colorbar(cax, ticks=[min(e), (min(e)+max(e))/2,max(e)])
cbar.ax.set_yticklabels([str(round(min(e),5)),
                         str(round((min(e)+max(e))/2,5)), round(max(e),5)])
# vertically xoriented colorbar
#plt.savefig(str(number)+".png")
#plt.show()