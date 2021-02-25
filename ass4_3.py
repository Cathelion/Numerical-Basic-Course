#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 16:49:06 2020

@author: flo
"""
import numpy as np
from scipy.integrate import quad 
import matplotlib.pyplot as plt

#make a polynomial class that is initialized with the 
#coefficients in form of a vector
# then the function(the polynomial) can be evaluated at every point
# just like a normal function definition (needed for quad too)
class polynomial:
    def __init__(self,v):
        self.v = v
        self.dim = len(v)-1
    
    # evaluates the funcntion by vector vector multiplication
    def evaluate(self,x):
        xDat = np.array([x**i for i in range(self.dim+1)])
        result = np.sum(np.multiply(self.v,xDat))
        return result       


# function to find the order of a method
def findOrder(method,a,b):
    order = -1   #just some starting value
    
    print("\n\n***now testing ",method.__name__, " ***" )
    
    # check if polynomials of increasing degree can be integrated
    for i in np.arange(1,10):
        print("\n---Test polynomials of degree ", i-1)
        flag = True
        
        # take 10 random polynomials of degree i 
        # so we avoid lucky cases
        for k in range(10):
            testcoeff = np.random.uniform(low=-1000, high =1000, size=i)
            poly = polynomial (testcoeff)
            num_int = method(poly.evaluate,a,b)
            num_real = quad(poly.evaluate,a,b)[0]
            error = np.abs(num_int-num_real )
            #print("error " , error )
            
            # if the error is to big - not equal anymore, break
            if (error >= 1e-10):
                flag = False
                order = i-1
                break
        
        if( not flag):
            break
        
    print(">>> order probably ", order)
    
    
# three different methods
def midpoint(f,a,b):
    return f((a+b)/2)*(b-a)

def leftpoint(f,a,b):
    return f(a)*(b-a)

def gauss2(f,a,b):  # a needs to be -1, b needs to be 1 
    # the c,b values are given by theoretical computation
    c_1 = -1/np.sqrt(3)
    c_2 = 1/np.sqrt(3)

    b_1 =1
    b_2 =1
    
    result = (b_1*f(c_1) + b_2*f(c_2))
    return result
    

# the three composite methods
def midpointComp(f,a,b,n):
    real_int = quad(f,a,b)[0]
    h = (b-a)/n
    result = 0
    for i in range(n):
        result = result + midpoint(f,a+i*h,a+(i+1)*h)
    return np.abs(result-real_int)

def leftpointComp(f,a,b,n):
    real_int = quad(f,a,b)[0]
    h = (b-a)/n
    result = 0
    for i in range(n):
        result = result + leftpoint(f,a+i*h,a+(i+1)*h)
    return np.abs(result-real_int)

def gauss2Comp(f,a,b,n):
    real_int = quad(f,a,b)[0]
    h = (b-a)/n
    result =0
    for j in range(n):
        x_j_1 = a+j*h
        x_j = a+(j+1)*h
        g = g_j(x_j_1,x_j,h,f)
        I_j = gauss2(g.evaluate,-1,1)
        result = result + I_j
    result = result /2 * h
    return np.abs(real_int-result)

class g_j:
    def __init__(self,x_j_1,x_j,h,f):
        self.x_j_1 = x_j_1
        self.x_j = x_j
        self.h =h
        self.f =f
        
    def evaluate(self,t):
        return self.f( 0.5*(self.x_j_1+self.x_j) + 0.5*self.h*t)


# some test functions
def f1(x):
    return np.sqrt(np.abs(x))

def f2(x):
    return np.exp(x**2)

def f3(x):
    if (x != 0):
        return x/np.abs(x)
    else:
        return 0

# a) part
findOrder(midpoint,-1,1)
findOrder(leftpoint,-1,1)
findOrder(gauss2,-1,1)


# c) part
xDat = list(range(1,20))

yDat11 = [midpointComp(f1,-1,1,n) for n in xDat]
yDat12 = [leftpointComp(f1,-1,1,n) for n in xDat]
yDat13 = [gauss2Comp(f1,-1,1,n) for n in xDat]

yDat21 = [midpointComp(f2,0,1,n) for n in xDat]
yDat22 = [leftpointComp(f2,0,1,n) for n in xDat]
yDat23 = [gauss2Comp(f2,0,1,n) for n in xDat]

yDat31 = [midpointComp(f3,-1,1,n) for n in xDat]
yDat32 = [leftpointComp(f3,-1,1,n) for n in xDat]
yDat33 = [gauss2Comp(f3,-1,1,n) for n in xDat]

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))
axes[0].plot(xDat, yDat11, label = "midpoint" ,color='r')
axes[0].plot(xDat, yDat12, label = "leftpoint" ,color='b')
axes[0].plot(xDat, yDat13, label = "gauss2",color='g')
axes[0].set_title('sqrt(|x|)')
axes[0].set_xlabel('number of subintervals')
axes[0].set_ylabel('abs error')
axes[0].legend()

axes[1].plot(xDat, yDat21, label = "midpoint",color='r')
axes[1].plot(xDat, yDat22, label = "leftpoint" ,color='b')
axes[1].plot(xDat, yDat23, label = "gauss2",color='g')
axes[1].set_title('exp(x^2)')
axes[1].set_xlabel('number of subintervals')
axes[1].set_ylabel('abs error')
axes[1].legend()

axes[2].plot(xDat, yDat31, label = "midpoint",color='r')
axes[2].plot(xDat, yDat32, label = "leftpoint" ,color='b')
axes[2].plot(xDat, yDat33, label = "gauss2" ,color='g')
axes[2].set_title('sign(x)')
axes[2].set_xlabel('number of subintervals')
axes[2].set_ylabel('abs error')
axes[2].legend()


print("Example for error spike; n=3")
h_ = 2/3
intervals = [-1+i*2/3 for i in range(4)]
print(intervals)

print("Example for error spike; n=9")
h_ = 2/9
intervals = [-1+i*2/9 for i in range(10)]
print(intervals)
