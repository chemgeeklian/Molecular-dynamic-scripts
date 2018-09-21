# -*- coding: utf-8 -*-
"""
Spyder Editor

By chemgeeklian, 2018/9/19
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import, error if not imported

class randwalk:
    def __init__(self,dim,norm,iteration):
        self.d=int(dim)
        self.r=norm
        self.i=iteration
    def direct(dimension,distance):
       rd=np.random.rand(dimension-1) 
       phai=rd*np.pi
       phai[-1]*=2
       
       direction=np.zeros(dimension)
       direction[0]=np.cos(phai[0]) #https://en.wikipedia.org/wiki/N-sphere
       sins=np.sin(phai[0])
       for i in range(1,dimension-1):
           direction[i]=sins*np.cos(phai[i])
           sins=sins*np.sin(phai[i-1])
       direction[-1]=sins
       return(distance*direction)

setpara=randwalk(3,1.0,500)
track=np.zeros([setpara.i,setpara.d]) #start from original point

for num in range(setpara.i):
    walk=randwalk.direct(setpara.d,setpara.r)
    track[num]=track[num-1]+walk
track=track.reshape([-1,setpara.d])

#specially to plot 3D track
fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(track[:,0],track[:,1],track[:,2], lw=1)
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
ax.set_title("3D Projection of Brownian Movement")

plt.show()
