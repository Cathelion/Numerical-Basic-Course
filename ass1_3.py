#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 21:03:35 2020

@author: flo
"""

import numpy as np
import matplotlib.pyplot as plt

class interpolate:
    
    def __init__(self,x_list, y_list):
        self.dim = -1
        self.dataX = x_list
        self.dataY = y_list
        self.polynomial = np.zeros(0)
        self.w = np.zeros(1)
        self.w[0] = 1
    
    # help function to multiply two polynomials if their product is in P_n
    def multiply_pol(self,v1,v2):
        res = np.zeros(len(v1))
        for i in np.arange(len(v2)):
            help_arr = np.roll(v1,i) * v2[i]   # shift 
            res = res+ help_arr
        return res
    
    # evaluates a polynomial in vector form at a given x
    def evaluate_pol(self,v1,x):
        res = 0
        for i in range(len(v1)):
            res = res + v1[i] * x**i
        return res
    
    # add each point of the given list
    def compute_poly(self):
        for i in range(len(self.dataX)):
            self.add_point(self.dataX[i], self.dataY[i])

        
    def add_point(self, newX, newY):
        self.dim = self.dim +1      # just keeping track of the dimension
        
        # compute factor c_n
        c_currC = (newY - self.evaluate_pol(self.polynomial, newX))/ self.evaluate_pol(self.w,newX)
        
        # help array so that mult. is possible
        helpArr = np.zeros(len(self.polynomial)+1)
        helpArr[:-1] = self.polynomial
        
        # compute new polynomial
        newP = helpArr + c_currC * self.w
        
        # make new w
        poly_factor = np.zeros(2)
        poly_factor[0] = - newX
        poly_factor[1] = 1
        helpArr2 = np.zeros(len(self.w)+1)
        helpArr2 [:-1] = self.w
        newW = self.multiply_pol(helpArr2 , poly_factor)
        
        # set field values
        self.w = newW
        self.polynomial = newP
    
    
    def print(self):
        return self.polynomial
        
    
    def plot(self):
        t_vals = np.arange(-2,3,0.01)
        
        y_vals = [self.evaluate_pol(self.polynomial, t) for t in t_vals]
        
        fig, ax = plt.subplots()
        ax.plot(t_vals, y_vals)
        ax.grid()
        plt.show()
        
        
        
x = [-1,0,2,3]
y=[2,6,4,30]
newton1 = interpolate(x,y)
newton1.compute_poly()
print(newton1.print())
newton1.add_point(1,5)
print(newton1.print())
newton1.plot()
newton1.add_point(4,8)
print(newton1.print())
newton1.plot()
