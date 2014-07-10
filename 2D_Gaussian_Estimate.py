# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 19:21:46 2014

@author: luis
"""

## This file draws from a bivariate gaussian distribution 
## and attempts to find the parameters drawn from using
## curve fit. Plots the resultant gaussian contours.

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import gauss
import plotutils

# Produce a number of points in x-y from 1 2D distribution. 
mean = [0,3]
cov = [[2,-20],[-20,2]]
N = 5000
x,y = np.random.multivariate_normal(mean,cov,N).T
plt.scatter(y,x,marker='x',linewidths=0.5) # Note transpose switches y and x

# Prep bins for histogram
bin_size = 0.2
max_edge = 4*(np.sqrt(cov[0][0])+np.sqrt(cov[1][1])) 
min_edge = -max_edge
bin_num = (max_edge-min_edge)/bin_size
bin_numPlus1 = bin_num + 1
bins = np.linspace(min_edge,max_edge,bin_numPlus1)

# Produce 2D histogram
H,xedges,yedges = np.histogram2d(x,y,bins,normed=False)
bin_centers_x = (xedges[:-1]+xedges[1:])/2.0
x_len = len(bin_centers_x)
bin_centers_y = (yedges[:-1]+yedges[1:])/2.0
y_len = len(bin_centers_y)
X,Y = np.meshgrid(bin_centers_x,bin_centers_y)

H = H.ravel()

# Initial Guess
p0 = (H.max(),mean[0],mean[1],cov[0][0],cov[1][1],0.5,0)

# Curve Fit parameters
coeff, var_matrix = curve_fit(gauss.mult_gaussFun_Fit,(X,Y),H,p0=p0)

Z = gauss.mult_gaussFun_Fit((X,Y),*coeff)
Z = Z.reshape((x_len,y_len))
plt.contour(X,Y,Z)
plt.colorbar()
plt.show()