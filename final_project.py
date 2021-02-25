#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 21:02:14 2020

@author: flo
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

h=0.1
g = 9.81

# help class
class onestep:
    def __init__(self, y1,y2):
        self.y1 =y1
        self.y2 =y2
        
    # define equation of y1 implicitly
    def eq (self,y1Curr):
        y1 = self.y1
        y2 = self.y2
        return y1 + h/2 * (y2+y2 - h*g/2*(np.sin(y1) + np.sin(y1Curr)))-y1Curr
    
    # solve y1 and then get y2 for free
    def solver(self):
        y1Curr = fsolve(self.eq,1)[0]
        y2Curr = (y1Curr-self.y1)*2/h - self.y2
        return y1Curr, y2Curr
    
    def reset(self, y1,y2):
        self.y1 = y1
        self.y2 = y2
    

# starting 
y1Curr, y2Curr = (np.pi/2 , 0)
y1L=[]
y2L=[]

# one instance of class needed
currantStep = onestep(y1Curr, y2Curr)

for i in range(int(5/h)):
    print(y1Curr,y2Curr)
    y1L.append(y1Curr)
    y2L.append(y2Curr)
    currantStep.reset(y1Curr, y2Curr)
    y1New,y2New = currantStep.solver()
    y1Curr = y1New
    y2Curr = y2New
    
    
plt.figure(figsize = (10,5))
plt.plot(y1L,y2L,label="phasespace")
x = np.linspace(0,5,int(5/h))
plt.plot(x,y1L, label="angle $alpha$")
plt.plot(x,y2L, label="velocity")
plt.legend()







