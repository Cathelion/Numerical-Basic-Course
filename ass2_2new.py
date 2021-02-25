#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 12:24:27 2020

@author: flo
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

def cubspline(xint, yint, xVal):   #xint, yint have length m+1
    
    # set constants
    m = len(xint)-1
    h = xint[1] - xint[0]
    
    # y-vector and 141 matrix
    y_vec =  6/h**2 * np.diff(np.diff((yint))) 
    A = np.diag(np.ones(m-2),-1) + np.diag((m-1)*[4],0) + np.diag(np.ones(m-2),1)
    
    # solve for sigma_1 until sigma_(m-1)
    sigma = linalg.solve(A,y_vec)
    
    # add sigma_0 and sigma_m (which are both zero from the natural boundary condition)
    tmp = np.zeros(len(sigma)+2)
    tmp[1:-1] = sigma
    sigma = tmp  #now sigma has length m+1
    
    # set up coeff matrix
    Mat = np.zeros([m,4])
    
    # compute coeff_i
    for i in np.arange(0,m):
        list_i = getCoeff(sigma[i] , sigma[i+1], yint[i],yint[i+1] ,h)
        Mat[i] =  list_i
    
    
    # find interval of xVal and thereby corresponding polynomial
    for j in np.arange(m):
        if (xint[j] <= xVal and xVal <=xint[j+1]):
            break
    
    
    spline = Mat[j]  # extract corresponding coefficients
    x_i = xint[j]
    yVal = spline[0]*(xVal - x_i)**3 + spline[1]*(xVal - x_i)**2 + spline[2]*(xVal - x_i) + spline[3]
    
    return yVal

def getCoeff(sigma_i, sigma_i1, y_i,y_i1,h):
    d = y_i
    b = sigma_i /2
    a= (sigma_i1 - sigma_i) / (6* h)
    c = (y_i1 - y_i)/h - h*(2*sigma_i + sigma_i1)/6
    
    return [a,b,c,d]


f1 = plt.figure(1, figsize = (10,10))

# f(t) = exp(-4t^2) , t in [-1,1] 
def f(t):
    return np.e**(-4*t**2)
xX1 = np.linspace(-1,1,200)
yY1 = [f(x) for x in xX1]

plt.subplot(211)
plt.plot(xX1,yY1, label="f(x)")
plt.grid()

# 5 nodes
xint_1 = np.linspace(-1,1,5)
yint_1 = [f(x) for x in xint_1]
yY_1cs = [cubspline(xint_1,yint_1,xVal) for xVal in xX1]
plt.plot(xX1,yY_1cs, label="5 nodes")

# 12 nodes
xint_2 = np.linspace(-1,1,12)
yint_2 = [f(x) for x in xint_2]
yY_2cs = [cubspline(xint_2,yint_2,xVal) for xVal in xX1]
plt.plot(xX1,yY_2cs, label="12 nodes")
plt.legend()

# g(t) = 1/(1+25t^2) , t in [-1,1] 
def g(t):
    return 1/(1+25*t**2) 
yY2 = [g(x) for x in xX1]

plt.subplot(212)
plt.plot(xX1,yY2, label="g(x)")
plt.grid()

# 15 nodes
xint_3 = np.linspace(-1,1,15)
yint_3 = [g(x) for x in xint_3]
yY_3cs = [cubspline(xint_3,yint_3,xVal) for xVal in xX1]
plt.plot(xX1,yY_3cs, label="15 nodes")

# 21 nodes
xint_4 = np.linspace(-1,1,21)
yint_4 = [g(x) for x in xint_4]
yY_4cs = [cubspline(xint_4,yint_4,xVal) for xVal in xX1]
plt.plot(xX1,yY_4cs, label="21 nodes")
plt.legend()


