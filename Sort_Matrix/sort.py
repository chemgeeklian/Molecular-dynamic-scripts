# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 15:16:39 2017

Automatically sort the coordination matrix of a molecular

@author: Xinran Lian
"""

import pandas as pd
import numpy as np

k = 0 
'''
k is used to set the origin point of the molecularv, 
it can be manually changed to avoid symmetry issues.
'''
data = pd.read_table('initial.txt',header=None,delim_whitespace=True,index_col=0)
a = np.array(data)
data['4'] = 0

for i in range(0,a.shape[0]):
    d = np.sqrt((a[i][0]-a[k][0])**2+(a[i][1]-a[k][1])**2+(a[i][2]-a[k][2])**2)
    data.iloc[i,3] = d

data1 = data.sort_values(by='4',ascending=True,kind='quicksort')
output = data1.iloc[:,0:3]
#Sort the atoms by distance between them and the zero point

f = open('out.txt','w')
for k in range(0,a.shape[0]):
    f.write(str(output.index[k]))
    f.write("   ")
    for m in range(0,3):
        f.write(str(output.iloc[k,m]))
        f.write("   ")
    f.write('\n')
f.close()
