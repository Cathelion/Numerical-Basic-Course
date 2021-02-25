#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 20:48:05 2020

@author: flo
"""

from scipy import linalg
import numpy as np
import matplotlib.pyplot as plt
from s1002 import s1002


""" exact same as in exercise 2"""
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

# helper function
def getCoeff(sigma_i, sigma_i1, y_i,y_i1,h):
    d = y_i
    b = sigma_i /2
    a= (sigma_i1 - sigma_i) / (6* h)
    c = (y_i1 - y_i)/h - h*(2*sigma_i + sigma_i1)/6
    
    return [a,b,c,d]


xDat = np.linspace(-70,60,200)
yDat = [-s1002(x) for x in xDat]

def getSpineDat(n):
    xD = np.linspace(-70,60,n)
    yD = [-s1002(x) for x in xD]
    return xD,yD

startVals5 = getSpineDat(5)
startVals15 = getSpineDat(15)
startVals50 = getSpineDat(50)
yDatC5 = [cubspline(startVals5[0], startVals5[1], x) for x in xDat]
yDatC15 = [cubspline(startVals15[0], startVals15[1], x) for x in xDat]
yDatC50 = [cubspline(startVals50[0], startVals50[1], x) for x in xDat]
plt.plot(xDat,yDat)
plt.plot(xDat,yDatC5,label="5 nodes")
plt.plot(xDat,yDatC15,label="15 nodes")
plt.plot(xDat,yDatC50,label="50 nodes")
plt.legend()
plt.grid()





