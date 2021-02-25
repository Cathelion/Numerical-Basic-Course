#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 21:59:04 2020

@author: flo
"""
import numpy as np
from scipy import linalg

def cubspline(xint, yint):   #xint, yint have length m+1
    
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
    
    return Mat

def getCoeff(sigma_i, sigma_i1, y_i,y_i1,h):
    d = y_i
    b = sigma_i /2
    a= (sigma_i1 - sigma_i) / (6* h)
    c = (y_i1 - y_i)/h - h*(2*sigma_i + sigma_i1)/6
    
    return [a,b,c,d]


# f(t) = exp(-4t^2) , t in [-1,1] 
def f(t):
    return np.e**(-4*t**2)

# 5 nodes
xint_1 = np.linspace(-1,1,5)
yint_1 = [f(x) for x in xint_1]
coeff_1 = cubspline(xint_1,yint_1)
print("with 5 nodes:" ,coeff_1)

# 12 nodes
xint_2 = np.linspace(-1,1,12)
yint_2 = [f(x) for x in xint_2]
coeff_2 = cubspline(xint_2,yint_2)
print("with 12 nodes:" ,coeff_2)


# g(t) = 1/(1+25t^2) , t in [-1,1] 
def g(t):
    return 1/(1+25*t**2) 

# 15 nodes
xint_3 = np.linspace(-1,1,15)
yint_3 = [g(x) for x in xint_3]
coeff_3 = cubspline(xint_3,yint_3)
print("with 15 nodes:" ,coeff_3)

# 21 nodes
xint_4 = np.linspace(-1,1,21)
yint_4 = [g(x) for x in xint_4]
coeff_4 = cubspline(xint_4,yint_4)
print("with 21 nodes:" ,coeff_4)








