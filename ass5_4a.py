#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 21:33:23 2020

@author: flo
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def y_1(t):
    return 100/99 * np.exp(-t) - 1/99 * np.exp(-100*t)

def y_2(t):
    return np.exp(-100*t)

def explicit(M,start,n,a,b):
    h = (b-a)/n
    vectors = []
    y_currant = start
    vectors.append(y_currant)
    for i in range(n):
        y_new = y_currant + h* np.dot(M,y_currant)
        vectors.append(y_new)
        y_currant = y_new
    return np.array(vectors).T

def implicit(M,start,n,a,b):
    h = (b-a)/n
    vectors = []
    y_currant = start
    vectors.append(y_currant)
    I = np.eye(2)
    A = I-h*M
    for i in range(n):
        y_new = np.linalg.solve(A,y_currant)
        vectors.append(y_new)
        y_currant = y_new
    return np.array(vectors).T

def err_ex(M,start,n):
    res = explicit(M,start,n,0,1)
    xDat = np.linspace(0,1,n+1)
    yDat00 = [y_1(t) for t in xDat]
    yDat11 = [y_2(t)for t in xDat]
    return np.sum(np.abs(res[0] - yDat00))  +  np.sum(np.abs(res[1] - yDat11))

def err_im(M,start,n):
    res = implicit(M,start,n,0,1)
    xDat = np.linspace(0,1,n+1)
    yDat00 = [y_1(t) for t in xDat]
    yDat11 = [y_2(t)for t in xDat]
    return np.sum(np.abs(res[0] - yDat00))  +  np.sum(np.abs(res[1] - yDat11))

M = np.array([[-1,1],[0,-100]])
start = np.array([1,1])
xDat = np.linspace(0,1,200)
yDat1 = [y_1(t) for t in xDat]
yDat2 = [y_2(t) for t in xDat]


plt.figure(figsize=(15,10))
plt.xlabel("y_1(t)")
plt.ylabel("y_2(t)")
"""
plt.plot(yDat1,yDat2,label ="real")

res1 = explicit(M, start, 70, 0, 1)
plt.plot(res1[0],res1[1], label="numeric 70")

res2 = explicit(M, start, 150, 0, 1)
plt.plot(res2[0],res2[1], label="numeric 150")

res2 = explicit(M, start, 250, 0, 1)
plt.plot(res2[0],res2[1], label="numeric 250")

plt.legend()
"""

plt.plot(yDat1,yDat2,label ="real y")
res2 = implicit(M, start, 12, 0, 1)
plt.plot(res2[0],res2[1], label="numeric 12")


res3 = implicit(M, start, 50, 0, 1)
plt.plot(res3[0],res3[1], label="numeric 50")

res4 = implicit(M, start, 100, 0, 1)
plt.plot(res4[0],res4[1], label="numeric 100")
plt.legend()


"""
nDat = [10*k for k in range(1,10)]
yDat1 = [err_ex(M,start,n) for n in nDat ]
yDat2 = [err_im(M,start,n) for n in nDat ]
#plt.plot(nDat,yDat1,label="ex")
plt.plot(nDat,yDat2,label="im")
plt.legend()
"""
"""
n =150
res1 =explicit(M,start,n,0,1)
res2 =implicit(M,start,n,0,1)

xDat = np.linspace(0,1,n+1)
xDat2 = np.linspace(0,1,200)
yDat0 = res1[0]
yDat1 = res1[1]
yDat2 = res2[0]
yDat3 = res2[1]
yDat00 = [y_1(t) for t in xDat2]
yDat11 = [y_2(t)for t in xDat2]

plt.figure(figsize = (15,10))
#plt.plot(xDat,yDat0,label="numeric y_1")
#plt.plot(xDat,yDat1,label="numeric y_2")
plt.xlabel("t")
plt.ylabel("y(t)")
plt.plot(xDat,yDat2,label=" numeric y_1")
plt.plot(xDat,yDat3,label="numeric y_2")
plt.plot(xDat2,yDat00,label="y_1")
plt.plot(xDat2,yDat11,label="y_2")
plt.legend()
plt.grid()
"""