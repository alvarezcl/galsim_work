# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 10:19:14 2014

@author: luis
"""

## This file loops through different sheared galsim objects
## in the form of gaussians and calculates the gaussian
## corresponding to the data using lmfit. Ellipticity
## estimates are then found and compared to the true value.

import numpy as np
import lmfit
import matplotlib.pyplot as plt
import gauss
import plotutils
import galsim
from mpl_toolkits.mplot3d import Axes3D

mag_shear = 0.5
# Ellipticities
el = [[mag_shear,0.0],[mag_shear,mag_shear],[0.0,mag_shear],[-mag_shear,mag_shear],
      [-mag_shear,0.0],[-mag_shear,-mag_shear],[0.0,-mag_shear],[mag_shear,-mag_shear]]
pi = np.pi

# Store [e1_est/e1_galsim,e2_est/e2_galsim]
stats = []
angle = np.linspace(0,7*np.pi/8,8)
count = 0
for val in el:

    # Parameters
    gal_flux = 1e5     # total counts on the image
    gal_HLR = 8.       # arcsec
    psf_sigma = 1.     # arcsec
    pixel_scale = 0.2  # arcsec / pixel
    noise = 0.        # standard deviation of the counts in each pixel
    size_image = (100,100) # Pixels or Arcsec?
    
    # Shear Param
    e1 = val[0]
    e2 = val[1]
    
    # Create primary galaxy 
    galaxy = galsim.Gaussian(flux=gal_flux,half_light_radius=gal_HLR)

    # Provide the shear
    galaxy = galaxy.shear(e1=e1,e2=e2)
    
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
    
#    plt.imshow(H,origin='lower')
#    plt.xlabel('$Pixels$'); plt.ylabel('$Pixels$')
#    plt.colorbar()
#    plt.title('Galaxy distribution witth ellipticity contours at ' + r'$\frac{%d\pi}{8}$'%(count),fontsize=18)
        
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
    
    p0 = (H.max(),x_cen,y_cen,30,30,0)
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
    
# Attempt to plot figures of snapshots    
#    if count == 0:
#        ax1 = fig.add_subplot(131)
#        a = ax1.imshow(H,interpolation='none',origin='lower')
#        plt.title('Angle: ' + r'$\frac{%d\pi}{8}$'%(count),fontsize=18)
#        plt.colorbar(a,shrink=0.4)
#        plt.xlabel('$Pixels$'); plt.ylabel('$Pixels$')
#        plt.contour(X,Y,mult_gaussFun_Fit((X,Y),*p_est))
#
#
#    if count == 5:
#        ax2 = fig.add_subplot(132)
#        a = ax2.imshow(H,interpolation='none',origin='lower')
#        plt.title('Angle: ' + r'$\frac{%d\pi}{8}$'%(count),fontsize=18)
#        plt.colorbar(a,shrink=0.4)
#        plt.xlabel('$Pixels$'); plt.ylabel('$Pixels$')
#        plt.contour(X,Y,mult_gaussFun_Fit((X,Y),*p_est))
#
#    
#    if count == 8:
#        ax3 = fig.add_subplot(131)
#        a = ax3.imshow(H,interpolation='none',origin='lower')
#        plt.title('Angle: ' + r'$\frac{%d\pi}{8}$'%(count),fontsize=18)
#        plt.colorbar(a,shrink=0.4)
#        plt.xlabel('$Pixels$'); plt.ylabel('$Pixels$')
#        plt.contour(X,Y,mult_gaussFun_Fit((X,Y),*p_est))
#        plt.show()
    
    count += 1    
    divisor = sigma_x_est**2 + sigma_y_est**2 + 2*(sigma_x_est*sigma_y_est - sigma_xy_est**2)**(1/2)
    e1_est = (sigma_x_est**2 - sigma_y_est**2)/divisor
    e2_est = 2*sigma_xy_est/divisor
    
#    results = im.FindAdaptiveMom()
#    shear = results.observed_shape
#    e1_galsim = shear.e1
#    e2_galsim = shear.e2
    
    stats.append([(e1_est-e1),(e2_est-e2)])
    
stats = np.array(stats)

plt.scatter(angle,stats[:,0],c='r',marker='x',label='$e_{1est}-e_{1true}$')
plt.scatter(angle,stats[:,1],label='$e_{2est}-e_{2true}$')
plt.xticks([0, pi/8, 2*pi/8, 3*pi/8, 4*pi/8, 5*pi/8, 6*pi/8, 7*pi/8],
           ['$0$', r'$\frac{\pi}{8}$', r'$\frac{2\pi}{8}$', r'$\frac{3\pi}{8}$', r'$\frac{4\pi}{8}$',
           r'$\frac{5\pi}{8}$', r'$\frac{6\pi}{8}$', r'$\frac{7\pi}{8}$'],fontsize = 16)
plt.ylabel('Difference: $e_1, e_2$',fontsize=16)           
plt.legend(prop={'size':16})
plt.title('Measured Minus True Ellipticity vs Angle of Orientation')
plt.show()