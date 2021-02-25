#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 20:19:28 2020

@author: flo
"""

import numpy as np
import scipy.integrate as inte
import matplotlib.pyplot as plt

def simpson1(f, a,b):
    return (b-a) /6 * (f(a) + 4* f ( (a+b)/2 ) + f(b) )

def simpson2(f,a,b):   # just a little bit faster than the simpson1 executed twice
    return (b-a) / 12 * (f(a) + 4* f( a+(b-a)/4) + 2 * f( (a+b)/2 ) + 4*f(a+ 3*(b-a)/4) +f(b) )

def f(x):
    return np.sqrt(np.abs(x))

def int_approx(f,a,b, err, t_0, t_e,intervals):
    delta = 15 * err * (t_e -t_0)
    if np.abs(simpson1(f,a,b)-simpson2(f,a,b)) > (b-a)*delta:
        intervals.append((a+b)/2)
        left = int_approx(f,a,(a+b)/2,err,t_0,t_e,intervals)
        right = int_approx(f,(a+b)/2,b,err,t_0,t_e,intervals)
        return left+right
    else:
        return simpson1(f,a,b)
    
        
a = -1
b = 1
error = 0.000005
intervals = [a,b]
I_head = int_approx(f,a,b, error,-1,1,intervals)
realI = inte.quad(f,a,b)[0]
print("real Integral (quad): " , realI)
print("approximated Integral(simpson): " , I_head)
print("real error(simpson's integral - quad's integral) : " , np.abs(I_head-realI))
print("estimated error : ", error*(len(intervals)+1))
print(len(intervals),  intervals)
#intervals2 = [x  for x in intervals if np.abs(x) < 0.008]
plt.eventplot(intervals,orientation='horizontal', colors='b', linelengths =0.05, lineoffsets = -0.1, label="intervals")
xD = np.linspace(a,b,201)
yD = [f(x) for x in xD]
plt.plot(xD,yD, label= "f(x)", color = 'r')
plt.legend(loc= 'center right')





