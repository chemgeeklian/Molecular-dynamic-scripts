# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 13:54:58 2018
simulation of thermal desorption spectrometry (TDS)

@author: Xinran Lian
"""

from scipy.integrate import ode
import numpy as np
from matplotlib import pyplot as plt

y0,t0=1,0
T0=1
t1=50
dt=0.02
Elr=40 #activation energy
Erl=20

def f1(t,y):
    k=[np.exp(-Elr/(T0+t)),np.exp(-Erl/(T0+t))]
    f=k[0]*0.1-(k[0]+k[1])*y
    return(f)

def f0(t,y):
    k=[np.exp(-Elr/(T0+t)),np.exp(-Erl/(T0+t))]
    f=-k[0]*y
    return(f)

def solveode(f):
    r = ode(f).set_integrator('vode', method='bdf')
    r.set_initial_value(y0, t0)
    time=[t0]
    yy=[y0]
    while r.successful() and r.t < t1:
        r.integrate(r.t+dt)
        time.append(r.t)
        yy.append(float(r.y[0]))
    yybig=np.array(yy[1:])
    yysmall=np.array(yy[:-1])
    difference=(yybig-yysmall)/dt
    return time,difference

time,difference=solveode(f0)
time,difference1=solveode(f1)
time=np.array(time)
plt.plot(300+10*time[:-1],-difference,'b',label='forward only')
plt.plot(300+10*time[:-1],-difference1,'r',label='forward and backward')
plt.legend()
plt.xlabel('T')
plt.ylabel('-d[$M_{n}O_{m}^+$]/d[t]')
plt.title('demo of TDS')
plt.show()