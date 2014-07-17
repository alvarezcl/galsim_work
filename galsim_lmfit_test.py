# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 20:16:33 2014

@author: luis
"""

## This file simply uses lmfit on a galsim object in 
## an analytic gaussian form in order to obtain
## parameters of the best fit.

from __future__ import division
import numpy as np
import lmfit
import matplotlib.pyplot as plt
import gauss
import plotutils
import galsim

# Parameters
gal_flux = 1e5     # total counts on the image
gal_HLR = 8.       # arcsec
psf_sigma = 1.     # arcsec
pixel_scale = 0.2  # arcsec / pixel
noise = 0.        # standard deviation of the counts in each pixel

# Create primary galaxy 
galaxy = galsim.Gaussian(flux=gal_flux,half_light_radius=gal_HLR)

# Plot information from created galaxy on contour
im = galaxy.drawShoot(scale=pixel_scale)
bounds = im.bounds
center = im.center()
x_cen = center.x
y_cen = center.y
xmin = bounds.getXMin(); xmax = bounds.getXMax()
ymin = bounds.getYMin(); ymax = bounds.getYMax()
x = np.linspace(xmin,xmax,xmax-xmin+1)
y = np.linspace(ymin,ymax,ymax-ymin+1)
X,Y = np.meshgrid(x,y)
H = im.array
plt.imshow(H,interpolation='none',origin='lower')
plt.title('Galsim distribution with Gaussian Contour Fit by LM_Fit\n')
plt.colorbar()

# Obtain the estimate of the gaussian from lmfit

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

#p0 = np.random.uniform(0,10,6)
p0 = np.array([H.max(),x_cen,y_cen,34,34,0])*2.0
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
plt.xlabel('$Pixels$'); plt.ylabel('$Pixels$')
plt.show()
divisor = sigma_x_est**2 + sigma_y_est**2 + 2*(sigma_x_est*sigma_y_est - sigma_xy_est**2)**(1/2)
e1 = (sigma_x_est**2 - sigma_y_est**2)/divisor
e2 = 2*sigma_xy_est/divisor

results = im.FindAdaptiveMom()