# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 2017

@author: Xinran Lian
"""
import matplotlib.pyplot as plt
import numpy as np

out=[0,0,0,0,0,0,0,0,0,0]
for i in range(0,10):
    read="d"+str(i)+".log"
    rd=open(read,'r')
    linen = 1
    for line in rd:
        if not line.find('Keep R1 ints in memory in canonical form') == -1:
            a=rd.readlines(1)
        linen +=1
    rd.close()
    str1 = a[0]
    b=float(str1[20:35])
    out[i]=round(b, 6)
print(out)

k=[0,0,0,0,0,0,0,0,0,0]
for j in range(0,10):
    k[j]=-(-1.0-0.05*float(j))
print(k)

plt.plot(k, out, '-', lw=2)

plt.xlabel('C-H distance(a.u.)')
plt.ylabel('E(a.u.)')
plt.grid(True)
plt.axes().set_aspect('equal', 'datalim')
plt.ylim((-40.2,-40.1))
plt.xlim((1,1.45))

plt.show()