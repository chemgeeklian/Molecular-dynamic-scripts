# -*- coding: utf-8 -*-
"""
Created on Tue May  8 15:24:13 2018

@author: Xinran Lian
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
loc='wavepacket/z/z/'
number=200

packet=list(np.zeros(number))
time=np.load(loc+'timearray.npy')
coordzip=np.load('coordinate.npz')
u=coordzip['u']*1.8897261339213 #r1 + r2
v=coordzip['v']*1.8897261339213 #r1 - r2 convert armstrong to atomic unit
scaleu=np.size(u)
scalev=np.size(v)

for i in range(number):
    packet[i]=np.load(loc+str(i*100+1)+'.npy')

time_plot=np.empty(number)
sum0=np.zeros(number)
sum1=np.zeros(number)
pkt=np.empty([scaleu*scalev,number,2]) 
#axis0:spacial region ; axis1:time region; axis2:[0:ground; 1:first excited]

for i in range(number):
    #packet[i]=np.abs(packet[i])
    packet[i] = packet[i] / np.linalg.norm(packet[i])
    pkt[:,i,0]=packet[i][:scaleu*scalev]
    pkt[:,i,1]=packet[i][scaleu*scalev:]
    sum0[i]=np.sum(pkt[:,i,0]*pkt[:,i,0].conjugate())
    sum1[i]=np.sum(pkt[:,i,1]*pkt[:,i,1].conjugate())
    time_plot[i]=time[i*5]

test_wp=pkt[:,1,1].reshape([-1,scaleu]).T
fig, ax = plt.subplots(figsize=(8, 8*scaleu/scalev))
ax.set_title('Wavefunction of BeH2, ground state, t=T/50\n')
ax.set_xlabel('r1-r2 (a.u.)')
ax.set_ylabel('r1+r2 (a.u.)')
cax = ax.imshow(test_wp,extent=[-1*v[scalev-1], v[scalev-1], u[scaleu-1], u[0]], 
                interpolation='nearest', cmap=cm.coolwarm)
plt.savefig('movieground.png')
plt.show()


plt.figure()

plt.plot(time_plot, sum0)#,marker='o')
plt.plot(time_plot, sum1)#,marker='o')
plt.title('Wavepacket')
plt.xlabel('time (a.u.)')
plt.ylabel('Propotion')
plt.legend(('Ground','Excited1'))

#plt.savefig('Scan.png')
plt.show()