# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 2017

@author: Xinran Lian
"""
import subprocess #library needed to call external programs
import os
os.environ["GAUSS_EXEDIR"] = "D:\Gaussian\g09C01window\G09W"

for i in range(0,10):
    filename="c"+str(i)+".gjf"
    f=open(filename,'w')    # r只读，w可写，a追加
    testnote= """# hf/6-31g

CH4 test

0 1
 C 0.00000002 0.00000000 0.00000000
 H 0.35665894 -1.00881114 0.00000000
 H 0.35667056 0.50439737 0.87365134
 H 0.35667056 0.50439737 -0.87365134
 H """
    f.write(testnote)
    j=-1.0+0.05*float(-i)
    f.write(str(j))
    f.write(""" 0.00000000 0.00000000
            
    """)
    f.close()
    
gaussianexe='D:\Gaussian\g09C01window\G09W\g09.exe' #name of the gaussian program
for k in range(0,10):
    filename="c"+str(k)+".gjf"
    outpfile="d"+str(k)+".log"
    inputfile=filename
    outputfile=outpfile
    print("Calling Gaussian: "+gaussianexe+" "+inputfile+" "+outputfile)
    subprocess.call([gaussianexe,inputfile,outputfile]) #calling gaussian09. The output is written to outputfile

print("Done.")