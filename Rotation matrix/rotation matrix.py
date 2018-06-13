# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 10:19:49 2017

@author: Xinran Lian
"""

import pandas as pd
import numpy as np

i=0;j=1
data = pd.read_table('import.txt',header=None,delim_whitespace=True,index_col=0)
a = np.array(data)
c0 = np.array(pd.read_table('obj.txt',header=None,delim_whitespace=True,index_col=0))
c = c0+a[i]-c0[i]  

'''
read the molecular coordinations from txt.
i and j are two endpoint atoms of the main axis. i is the coordination of the 
atom which is fixed as the starting point (set as (0,0,0) here)
'''
b1 = np.array([a[j][0]-a[i][0],a[j][1]-a[i][1],a[j][2]-a[i][2]])
b2 = np.array([c[j][0]-a[i][0],c[j][1]-a[i][1],c[j][2]-a[i][2]]) 

'''
Create vectors of axis of the initial and rotated molecular, respectively.
And convert them to array for dot operation
The dot operation is for length of vectors, which is used for theta.
'''
theta = np.arccos(b1.dot(b2)/(np.sqrt(b1.dot(b1))*np.sqrt(b2.dot(b2))))
if theta == 0:
    print("No rotation necessary, exit(0)")
    exit(0)
elif theta == np.pi/2:
    axis = np.array([1,2,3])
else:
    axis = np.cross(b1,b2) #Normal vector of the plane determined by b1 and b2

lr = np.sqrt(axis[0]**2 + axis[1]**2) #length of shade of axis on xoy plane
if lr !=0:
    trans = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[-a[i][0],-a[i][1],-a[i][2],1]])
    r_yox = np.array([[axis[1]/lr,axis[0]/lr,0,0],[-axis[0]/lr,axis[1]/lr,0,0],[0,0,1,0],[0,0,0,1]])
    r_z = np.array([[1,0,0,0],[0,axis[2],lr,0],[0,-lr,axis[2],0],[0,0,0,1]])
    rot = np.array([[np.cos(theta), np.sin(theta),0,0],\
                   [-np.sin(theta),np.cos(theta),0,0],[0,0,1,0],[0,0,0,1]])
    rotmatrix = trans.dot(r_yox).dot(r_z).dot(rot).dot(np.linalg.inv(r_z).\
                  dot(np.linalg.inv(r_yox)).dot(np.linalg.inv(trans)))
else:
    rotmatrix = trans.dot(rot).dot(np.linalg.inv(trans))
    
print(rotmatrix)
'''
trans is the transition matrix moving the origin point to (0,0,0)
r_yox rotates the axis to yoz plane
r_z rotates the axis to make it coincide with z axis
    The net matrix make the molecular rotate around z axis, then inverse matrixes 
move it back to its original position.
Ref. http://blog.csdn.net/qweewqpkn/article/details/73142656 (Chinese)
'''

rst=np.zeros((a.shape[0],4),dtype=float) #creat output matrix
c = np.column_stack((c,np.ones(a.shape[0]).reshape(a.shape[0],1)))
print(c)
for k in range(0,a.shape[0]):
    rst[k] = np.dot(c[k],rotmatrix)#the matrix is input row by row

f = open('output.txt','w')
for k in range(0,a.shape[0]):
    f.write(str(data.index[k]))
    f.write("   ")
    for m in range(0,3):
        f.write(str(rst[k][m]))
        f.write("   ")
    f.write('\n')
f.close()