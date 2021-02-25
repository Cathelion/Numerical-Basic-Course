# -*- coding: utf-8 -*-
"""
Created on Thu May  7 16:16:09 2020

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt


def explicit(M,y_start,n,a,b):
    h = (b-a)/n
    vectors = []
    y_currant = y_start
    A_currant = np.array([[M[0][0] , M[0][1]*y_currant[0]] ,
                          [M[1][0] * y_currant[1] , M[1][1]]])
    vectors.append(y_currant)
    for i in range(n):
        
        # calculate y_new  formula of old values
        y_new = y_currant + h* np.dot(A_currant,y_currant)
        vectors.append(y_new)
        
        # make new matrix
        A_new = np.array([[M[0][0] , M[0][1]*y_new[0]] ,
                          [M[1][0] * y_new[1] , M[1][1]]])
        
        # update
        y_currant = y_new
        #print(y_currant)
        A_currant = A_new
    return np.array(vectors).T

days = 30
accuracy =100
start = np.array([1,1])
M = np.array([[1, -1.5],[1 , -10]])
res = explicit(M,start,days*accuracy,0,days)

plt.figure(figsize=(15,10))
plt.xlabel("x(t)")
plt.ylabel("y(t)")
plt.plot(res[0], res[1])
print (res[0][0], res[1][0])

