# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 14:50:00 2014

@author: luis
"""

# Based on drawLibrary.py  by LCA
# Small modifications by MSSG
# Start: 8/2/2014

## Library containing sets of functions for drawing galaxies and residuals
## for usage with fitters.

import galsim
import numpy as np
import gauss
import sys

# Draw from a distribution and return a binned image object with one galaxy 
# Shoots photons, see last line before it returns image -mg
def drawShoot_galaxy(flux, hlr, e1, e2, x0, y0, x_len, y_len, scale, func, seed):
    big_fft_params = galsim.GSParams(maximum_fft_size=10024000)    

    print " *********************** Using MSSG version of drawLibrary"

    if func is galsim.Gaussian:    
        gal = func(half_light_radius=hlr, flux=flux, gsparams=big_fft_params)

    if func is galsim.Sersic:                       
        deVauc_ix = 4 # deVauc bulge
        expl_ix = 1   # expl disk
        btodsize = 0.5
        gal_HLR = hlr
        gal_flux = flux
        gal  = galsim.Sersic(n=deVauc_ix, flux=gal_flux,half_light_radius=gal_HLR*btodsize) # bulge
        gal += galsim.Sersic(n=expl_ix, flux=gal_flux,half_light_radius=gal_HLR)  # add in expl

    # Same for either
    gal = gal.shear(g1=e1, g2=e2)
    gal = gal.shift(x0,y0)
    image = galsim.ImageD(x_len, y_len, scale=scale)
    image = gal.drawImage(image=image, method='phot',rng=seed)

    return image

# Use the analytic definition of an image profile for one galaxy 
# Uses FFT,  see last line before it returns image and compare to drawshoot_galaxy above  -mg
def draw_galaxy_1(flux, hlr, e1, e2, x0, y0, x_len, y_len, scale, galtype):
    big_fft_params = galsim.GSParams(maximum_fft_size=10024000)
    if galtype is galsim.Gaussian:        
        gal = galtype(half_light_radius=hlr, flux=flux, gsparams=big_fft_params)

    if galtype is galsim.Sersic:        
        deVauc_ix = 4 # deVauc bulge
        expl_ix = 1   # expl bulge
        btodsize = 0.5
        gal_HLR = hlr
        gal_flux = flux
        gal = galsim.Sersic(n=deVauc_ix, flux=gal_flux, half_light_radius=gal_HLR*btodsize) # bulge
        gal += galsim.Sersic(n=expl_ix, flux=gal_flux, half_light_radius=gal_HLR)  # add in expl

    # Same for either
    gal = gal.shear(g1=e1, g2=e2)
    gal = gal.shift(x0,y0)
    image = galsim.ImageD(x_len, y_len, scale=scale)
    image = gal.drawImage(image=image)
    
    return image

# The difference of the data and the model for one galaxy
# that is to be reduced.
def resid_1obj(param, target_image, x_len, y_len, scale, galtype):
    flux = param['flux_a'].value
    hlr = param['hlr_a'].value
    e1 = param['e1_a'].value
    e2 = param['e2_a'].value
    x0 = param['x0_a'].value
    y0 = param['y0_a'].value
    image = draw_galaxy_1(flux,hlr,e1,e2,x0,y0,x_len,y_len,scale, galtype)
    return (image-target_image).array.ravel()

# Draw from two distributions and return a binned image object with two 
# galaxies 
def drawShoot_2galaxies(flux_1,hlr_1,e1_1,e2_1,x_center1,y_center1,
                       flux_2,hlr_2,e1_2,e2_2,x_center2,y_center2,
                       x_len,y_len,scale,func_1,func_2,seed):
    big_fft_params = galsim.GSParams(maximum_fft_size=10024000)                       
    if func_1 is galsim.Gaussian:                           
        gal_1 = func_1(half_light_radius=hlr_1, flux=flux_1, gsparams=big_fft_params)
        gal_1 = gal_1.shear(g1=e1_1, g2=e2_1)
        gal_1 = gal_1.shift(x_center1,y_center1)
        image_1 = galsim.ImageD(x_len, y_len, scale=scale)
        image_1 = gal_1.drawImage(image=image_1,method='phot',rng=seed)
    if func_2 is galsim.Gaussian:    
        gal_2 = func_2(half_light_radius=hlr_2, flux=flux_2, gsparams=big_fft_params)
        gal_2 = gal_2.shear(g1=e1_2, g2=e2_2)
        gal_2 = gal_2.shift(x_center2,y_center2)
        image_2 = galsim.ImageD(x_len, y_len, scale=scale)
        image_2 = gal_2.drawImage(image=image_2,method='phot',rng=seed)        
    image = image_1 + image_2    
    return image

# Use the analytic definition of an image profile for two galaxies -- Gaussian or Sersic -mg
def draw_2galaxies(flux_1,hlr_1,e1_1,e2_1,x_center1,y_center1,
                  flux_2,hlr_2,e1_2,e2_2,x_center2,y_center2,
                  x_len,y_len,scale,func_1,func_2):
    big_fft_params = galsim.GSParams(maximum_fft_size=10024000)
    if func_1 is galsim.Gaussian:                       
        gal_1 = func_1(half_light_radius=hlr_1, flux=flux_1, gsparams=big_fft_params)
        gal_1 = gal_1.shear(g1=e1_1, g2=e2_1)
        gal_1 = gal_1.shift(x_center1,y_center1)

    if func_1 is galsim.Sersic:        
        deVauc_ix = 4 # deVauc bulge
        expl_ix = 1   # expl bulge
        btodsize = 0.5
        gal_HLR = hlr_1
        gal_flux = flux_1
        gal_1 = galsim.Sersic(n=deVauc_ix, flux=gal_flux, half_light_radius=gal_HLR*btodsize) # bulge
        gal_1 += galsim.Sersic(n=expl_ix, flux=gal_flux, half_light_radius=gal_HLR)  # add in expl

    image_1 = galsim.ImageD(x_len, y_len, scale=scale)
    image_1 = gal_1.drawImage(image=image_1)

    if func_2 is galsim.Gaussian:
        gal_2 = func_2(half_light_radius=hlr_2, flux=flux_2, gsparams=big_fft_params)
        gal_2 = gal_2.shear(g1=e1_2, g2=e2_2)
        gal_2 = gal_2.shift(x_center2,y_center2)


    if func_1 is galsim.Sersic:        
        deVauc_ix = 4 # deVauc bulge
        expl_ix = 1   # expl bulge
        btodsize = 0.5
        gal_HLR = hlr_2
        gal_flux = flux_2
        gal_2 = galsim.Sersic(n=deVauc_ix, flux=gal_flux, half_light_radius=gal_HLR*btodsize) # bulge
        gal_2 += galsim.Sersic(n=expl_ix, flux=gal_flux, half_light_radius=gal_HLR)  # add in expl

    image_2 = galsim.ImageD(x_len, y_len, scale=scale)
    image_2 = gal_2.drawImage(image=image_2)

    image = image_1 + image_2
    
    return image

# The difference of the data and the model for a pair of galaxies
# that is to be minimized by fitter - mg
def resid_2obj(param, target_image, x_len, y_len, scale, func_1, func_2):
    flux_1 = param['flux_1'].value
    hlr_1 = param['hlr_1'].value
    e1_1 = param['e1_1'].value
    e2_1 = param['e2_1'].value
    x_center1 = param['x_center1'].value
    y_center1 = param['y_center1'].value

    flux_2 = param['flux_2'].value
    hlr_2 = param['hlr_2'].value
    e1_2 = param['e1_2'].value
    e2_2 = param['e2_2'].value
    x_center2 = param['x_center2'].value
    y_center2 = param['y_center2'].value
    
    image1 = draw_galaxy_1(flux_1,hlr_1,e1_1,e2_1,x_center1,y_center1,x_len,y_len,scale,func_1)
    image2 = draw_galaxy_1(flux_2,hlr_2,e1_2,e2_2,x_center2,y_center2,x_len,y_len,scale,func_2)
    image = image1 + image2

    # Put just sum of one galaxy instead
    
    # Create error array
    error = np.sqrt(target_image.array.ravel())
    # Set the errors equal to 1 where 0 is.
    error[error==0] = 1
    return (image-target_image).array.ravel()
    
    

# Return the bounds of an image.
def return_domain(image):
    bounds = image.bounds
    center = image.center()
    x_cen = center.x
    y_cen = center.y
    xmin = bounds.getXMin(); xmax = bounds.getXMax()
    ymin = bounds.getYMin(); ymax = bounds.getYMax()
    x = np.linspace(xmin,xmax,xmax-xmin+1)
    y = np.linspace(ymin,ymax,ymax-ymin+1)
    X,Y = np.meshgrid(x,y)
    return x_cen,y_cen,x,y,X,Y
    
# Residual function to minimize for 2D gaussians.
def resid_gauss(params, data,(X,Y)):
    A = params['amplitude'].value
    x0 = params['x_mean'].value
    y0 = params['y_mean'].value
    sigma_x = params['sigma_x'].value
    sigma_y = params['sigma_y'].value
    sigma_xy = params['sigma_xy'].value
    model = gauss.mult_gaussFun_Fit((X,Y),*(A,x0,y0,sigma_x,sigma_y,sigma_xy))
    return (model - data).ravel()

# For plotting the contours of a circular profile
def contour_positive_circle(x,x0,y0,R):
    return y0 + np.sqrt(R**2 - (x-x0)**2)
    
# For plotting the contours of a circular profile
def contour_negative_circle(x,x0,y0,R):
    return y0 - np.sqrt(R**2 - (x-x0)**2)
    
