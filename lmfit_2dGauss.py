# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 15:50:34 2014

@author: luis
"""

## This file draws from a bivariate gaussian distribution 
## and attempts to find the parameters drawn from using
## lm fit. Plots the resultant gaussian contours on the
## drawn points and plots 2D histogram.

import numpy as np
import lmfit
import matplotlib.pyplot as plt
import gauss
import plotutils

# Produce a number of points in x-y from a 2D gaussian distribution.
# Note the transpose 
N = 1000
mean = [1,1]
cov = [[4,-1],[-1,3]]
sigma_x = np.sqrt(cov[0][0])
sigma_y = np.sqrt(cov[1][1])
sigma_xy = cov[0][1]
y,x = np.random.multivariate_normal(mean,cov,N).T

# Plot the scatter of the second draw
p2 = plt.scatter(y,x,marker='x',linewidths=0.5) # Note transpose switches y and x

# Prep bins for histogram
bin_size = 0.1
max_edge = 2*(np.sqrt(cov[0][0])+np.sqrt(cov[1][1])) 
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

# Return a gaussian distribution at an angle alpha from the x-axis
# from astroML for use with curve_fit
def mult_gaussFun_Fit((X,Y),*m):
    A,x0,y0,sigma_x,sigma_y,sigma_xy = m
    rho = sigma_xy/(sigma_x*sigma_y)
    a = 1/(2*(1-rho**2))
    z_sq = ((X-x0)/(sigma_x))**2 + ((Y-y0)/(sigma_y))**2 - 2*(sigma_xy/(sigma_x*sigma_y)**2)*(X-x0)*(Y-y0)
    Z = A*np.exp(-a*z_sq)
    return Z

# Residual function to minimize.
def resid(params, data,(X,Y)):
    A = params['amplitude'].value
    x0 = params['x_mean'].value
    y0 = params['y_mean'].value
    sigma_x = params['sigma_x'].value
    sigma_y = params['sigma_y'].value
    sigma_xy = params['sigma_xy'].value
    model = mult_gaussFun_Fit((X,Y),*(A,x0,y0,sigma_x,sigma_y,sigma_xy))
    return (model - data).ravel()

# Seed and Parameters

#p0 = (H.max(),mean[0],mean[0],sigma_x,sigma_y,sigma_xy)
p0 = np.random.uniform(0,10,6)
params = lmfit.Parameters()
params.add('amplitude',value=p0[0],min=0)
params.add('x_mean',value=p0[1])
params.add('y_mean',value=p0[2])
params.add('sigma_x',value=p0[3],min=0)
params.add('sigma_y',value=p0[4],min=0)
params.add('sigma_xy',value=p0[5])

# Extract the best-fit parameters
result = lmfit.minimize(resid,params,args=(H,(X,Y)))
lmfit.report_errors(result.params)

p_est = (result.params['amplitude'].value,result.params['x_mean'].value,result.params['y_mean'].value,
         result.params['sigma_x'].value,result.params['sigma_y'].value,result.params['sigma_xy'].value)
A_est = result.params['amplitude'].value
x0_est = result.params['x_mean'].value
y0_est = result.params['y_mean'].value
sigma_x_est = result.params['sigma_x'].value
sigma_y_est = result.params['sigma_y'].value
sigma_xy_est = result.params['sigma_xy'].value
plt.contour(X,Y,mult_gaussFun_Fit((X,Y),*p_est))
plt.legend([p2],['Points Drawn'])
plt.title('Scattered Points and Gaussian Contour Fit Using LM_Fit')
plt.xlabel('$x$'); plt.ylabel('$y$')
plt.show()

# Plot the 2d histogram of the points
fig = plotutils.hist_3dBar_points(y,x,bins)
ax = fig.gca(projection='3d')
ax.set_xlabel('$x$'); ax.set_ylabel('$y$'); ax.set_zlabel('$Frequency$')
ax.text2D(0.05,0.95,'3D Histogram of Point Density',transform=ax.transAxes)
plt.show()