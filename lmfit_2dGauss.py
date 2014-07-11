# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 15:50:34 2014

@author: luis
"""

import numpy as np
import lmfit
import matplotlib.pyplot as plt

# Produce a number of points in x-y from a 2D gaussian distribution.
# Note the transpose 
N = 1000
mean = [0,0]
cov = [[3,0],[0,1]]
x,y = np.random.multivariate_normal(mean,cov,N).T
# Plot the scatter of the second draw
plt.scatter(x,y,marker='x',linewidths=0.5) # Note transpose switches y and x

# Return a gaussian distribution at an angle alpha from the x-axis
# from astroML for use with curve_fit
def mult_gaussFun_Fit((X,Y),*m):
    A,x0,y0,sigma_x,sigma_y,sigma_xy = m
    rho = sigma_xy/(sigma_x*sigma_y)
    a = 1/(2*(1-rho**2))
    z_sq = ((X-x0)/(sigma_x))**2 + ((Y-y0)/(sigma_y))**2 - 2*(sigma_xy/(sigma_x*sigma_y)**2)*(X-x0)*(Y-y0)
    Z = A*np.exp(-a*z_sq)
    return Z
        
def resid(params, data,(X,Y)):
    A = params['amplitude'].value
    x0 = params['x_mean'].value
    y0 = params['y_mean'].value
    sigma_x = params['sigma_x'].value
    sigma_y = params['sigma_y'].value
    sigma_xy = params['sigma_xy'].value
    model = mult_gaussFun_Fit((X,Y),*(A,x0,y0,sigma_x,sigma_y,sigma_xy))
    return (model - data).ravel()

# Seed
p0 = (30,0,0,1,2,0.5)
params = lmfit.Parameters()
params.add('amplitude',value=p0[0])
params.add('x_mean',value=p0[1])
params.add('y_mean',value=p0[2])
params.add('sigma_x',value=p0[3])
params.add('sigma_y',value=p0[4])
params.add('sigma_xy',value=p0[5])


x = y = np.linspace(-15,15,1000)
truth_val = (15,1,1,2,3,0.25)
X,Y = np.meshgrid(x,y)
data = mult_gaussFun_Fit((X,Y),*truth_val)
data += np.random.randn(*data.shape)*0.1

result = lmfit.minimize(resid,params,args=(data,(X,Y)))