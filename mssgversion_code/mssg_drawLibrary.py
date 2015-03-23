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
import mssg_gauss
import sys
import ipdb
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

# Version of 10/14/2014
# Draw from a distribution and return a binned image object with one galaxy 
def draw_1Comp_galaxy(flux, hlr, e1, e2, x0, y0, x_len, y_len, scale, func, seed, fftORphotons = 'fft'):
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

    # Same for either Sersic or Gaussian
    gal = gal.shear(g1=e1, g2=e2)
    gal = gal.shift(x0,y0)
    image = galsim.ImageD(x_len, y_len, scale=scale)
    
    if fftORphotons == 'fft':
        image = gal.drawImage(image=image)
    else:
        image = gal.drawImage(image=image, method='phot',rng=seed)
    
    # Send it back    
    return image

# Version of 10/14/2014
# Draw from a distribution and return a binned image object with one galaxy 
def draw_2Comp_galaxy(flux, hlr, e1, e2, x0, y0, x_len, y_len, scale, seed, btodsize = 1, galtype = 'galsim.Gaussian', fftORphotons = 'fft', returnObjOrImg = 'img' ): # Default bulge radius = half disk radius
    big_fft_params = galsim.GSParams(maximum_fft_size=10024000)    

    print " *********************** Using MSSG version of drawLibrary"
    print " ** Making gal of type = " , galtype

    if galtype == 'galsim.Gaussian':    
        gal = galsim.Gaussian(half_light_radius=hlr, flux=flux, gsparams=big_fft_params)

    if galtype == 'galsim.Sersic':                       
        deVauc_ix = 4 # deVauc bulge
        expl_ix = 1   # expl disk

        gal_HLR = hlr
        gal_flux = flux
        gal  = galsim.Sersic(n=deVauc_ix, flux=gal_flux,half_light_radius=gal_HLR*btodsize) # bulge
        gal += galsim.Sersic(n=expl_ix, flux=gal_flux,half_light_radius=gal_HLR)  # add in expl

    # Same for either Sersic or Gaussian
    gal = gal.shear(g1=e1, g2=e2)
    gal = gal.shift(x0,y0)


    print " returnObjOrImg = " , returnObjOrImg 

    if returnObjOrImg == 'obj':
        return gal
    else:
        image = galsim.ImageD(x_len, y_len, scale=scale)
        if fftORphotons == 'fft':
            image = gal.drawImage(image=image)
        else:
            image = gal.drawImage(image=image, method='phot',rng=seed)

            # Send it back    
        return image



########################################## Older

# Draw from a distribution and return a binned image object with one galaxy 
# Shoots photons, see last line before it returns image -- that's only dif from FFT version -mg
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

    # Same for either Sersic or Gaussian
    gal = gal.shear(g1=e1, g2=e2)
    gal = gal.shift(x0,y0)
    image = galsim.ImageD(x_len, y_len, scale=scale)
    image = gal.drawImage(image=image, method='phot',rng=seed)

    return image

# Draw from a distribution and return a binned image object with one galaxy 
# Shoots photons, see last line before it returns image -mg
def drawShoot_2comp_galaxy(flux, hlr, e1, e2, x0, y0, x_len, y_len, scale, func, seed, btodsize = 0.5): # Default bulge radius = half disk radius
    big_fft_params = galsim.GSParams(maximum_fft_size=10024000)    

    print " *********************** Using MSSG version of drawLibrary"

    if func is galsim.Gaussian:    
        gal = func(half_light_radius=hlr, flux=flux, gsparams=big_fft_params)

    if func is galsim.Sersic:                       
        deVauc_ix = 4 # deVauc bulge
        expl_ix = 1   # expl disk

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
def draw_galaxy_1(flux, hlr, e1, e2, x0, y0, x_len, y_len, scale, galtype, dopsfconvln = 'n'):
    print " e1, e2 = ", e1,e2
    print " flux, hlr = ", flux, hlr 
    if e1 < -0.95:
        ipdb.set_trace()
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

    # Same for either Sersic or Gaussian
    gal = gal.shear(g1=e1, g2=e2)
    gal = gal.shift(x0,y0)

    if dopsfconvln == 'y':
        psfshr = 0.00
        psf = galsim.Moffat(beta=3, fwhm=0.85).shear(e1=psfshr,  e2 = -psfshr)

        gal = galsim.Convolve([gal,psf])

    protoimage = galsim.ImageD(x_len, y_len, scale=scale)
    image = gal.drawImage(image=protoimage)
    
    return image

# The difference of the data and the model for one galaxy
# that is to be reduced.
def resid_1obj(param, target_image, x_len, y_len, scale, galtype, dopsfconvln = 'n'):
    flux = param['flux_a'].value
    hlr = param['hlr_a'].value
    e1 = param['e1_a'].value
    e2 = param['e2_a'].value
    x0 = param['x0_a'].value
    y0 = param['y0_a'].value

#    ipdb.set_trace()

#    print " e1_a = ", e1_a 

    image = draw_galaxy_1(flux,hlr,e1,e2,x0,y0,x_len,y_len,scale, galtype, dopsfconvln)
    return (image-target_image).array.ravel()

# Draw from two distributions and return a binned image object with two 
# galaxies 
def drawShoot_2galaxies(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,
                       flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,
                       x_len,y_len,scale,func_a,func_b,seed):
    big_fft_params = galsim.GSParams(maximum_fft_size=10024000)                       
    if func_a is galsim.Gaussian:                           
        gal_a = func_a(half_light_radius=hlr_a, flux=flux_a, gsparams=big_fft_params)
        gal_a = gal_a.shear(g1=e1_a, g2=e2_a)
        gal_a = gal_a.shift(x0_a,y0_a)
        image_a = galsim.ImageD(x_len, y_len, scale=scale)
        image_a = gal_a.drawImage(image=image_a,method='phot',rng=seed)
    if func_b is galsim.Gaussian:    
        gal_b = func_b(half_light_radius=hlr_b, flux=flux_b, gsparams=big_fft_params)
        gal_b = gal_b.shear(g1=e1_b, g2=e2_b)
        gal_b = gal_b.shift(x0_b,y0_b)
        image_b = galsim.ImageD(x_len, y_len, scale=scale)
        image_b = gal_b.drawImage(image=image_b,method='phot',rng=seed)        
    image = image_a + image_b    
    return image

# Use the analytic definition of an image profile for two galaxies -- Gaussian or Sersic -mg
def draw_2galaxies(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,
                  flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,
                  x_len,y_len,scale,func_a,func_b):
    big_fft_params = galsim.GSParams(maximum_fft_size=10024000)
    if func_a is galsim.Gaussian:                       
        gal_a = func_a(half_light_radius=hlr_a, flux=flux_a, gsparams=big_fft_params)
        gal_a = gal_a.shear(g1=e1_a, g2=e2_a)
        gal_a = gal_a.shift(x0_a,y0_a)

    if func_a is galsim.Sersic:        
        deVauc_ix = 4 # deVauc bulge
        expl_ix = 1   # expl bulge
        btodsize = 0.5
        gal_HLR = hlr_a
        gal_flux = flux_a
        gal_a = galsim.Sersic(n=deVauc_ix, flux=gal_flux, half_light_radius=gal_HLR*btodsize) # bulge
        gal_a += galsim.Sersic(n=expl_ix, flux=gal_flux, half_light_radius=gal_HLR)  # add in expl

    image_a = galsim.ImageD(x_len, y_len, scale=scale)
    image_a = gal_a.drawImage(image=image_a)

    if func_b is galsim.Gaussian:
        gal_b = func_b(half_light_radius=hlr_b, flux=flux_b, gsparams=big_fft_params)
        gal_b = gal_b.shear(g1=e1_b, g2=e2_b)
        gal_b = gal_b.shift(x0_b,y0_b)


    if func_a is galsim.Sersic:        
        deVauc_ix = 4 # deVauc bulge
        expl_ix = 1   # expl bulge
        btodsize = 0.5
        gal_HLR = hlr_b
        gal_flux = flux_b
        gal_b = galsim.Sersic(n=deVauc_ix, flux=gal_flux, half_light_radius=gal_HLR*btodsize) # bulge
        gal_b += galsim.Sersic(n=expl_ix, flux=gal_flux, half_light_radius=gal_HLR)  # add in expl

    image_b = galsim.ImageD(x_len, y_len, scale=scale)
    image_b = gal_b.drawImage(image=image_b)

    image = image_a + image_b
    
    return image

# The difference of the data and the model for a pair of galaxies
# that is to be minimized by fitter - mg
def resid_2obj(param, target_image, x_len, y_len, scale, func_a, func_b):
#    ipdb.set_trace()
    flux_a = param['flux_a'].value
    hlr_a = param['hlr_a'].value
    e1_a = param['e1_a'].value
    e2_a = param['e2_a'].value
    x0_a = param['x0_a'].value
    y0_a = param['y0_a'].value

    flux_b = param['flux_b'].value
    hlr_b = param['hlr_b'].value
    e1_b = param['e1_b'].value
    e2_b = param['e2_b'].value
    x0_b = param['x0_b'].value
    y0_b = param['y0_b'].value

    dopsfconvln = 'y'
    
    image1 = draw_galaxy_1(flux_a, hlr_a, e1_a, e2_a, x0_a, y0_a, x_len, y_len, scale, func_a, dopsfconvln)
    image2 = draw_galaxy_1(flux_b, hlr_b, e1_b, e2_b, x0_b, y0_b, x_len, y_len, scale, func_b, dopsfconvln)
    image = image1 + image2

    # Put just sum of one galaxy instead
    
    # Create error array
#    error = np.sqrt(target_image.array.ravel())
    # Set the errors equal to 1 where 0 is.
#    error[error==0] = 1

    '''
    # Initialize the geometry of the grid
    gs = gridspec.GridSpec(1,3)
    fig = plt.figure(figsize=(10,5))

    # Plot the target data img
    f1 = fig.add_subplot(gs[0,0])
    plt.title('obj a')
    fig1 = f1.imshow(image1.array,interpolation='none',origin='lower')
    plt.colorbar(fig1, shrink=0.5) 

    # Plot the best fit img
    f2 = fig.add_subplot(gs[0,1])
    plt.title('obj b')
    fig2 = f2.imshow(image2.array,interpolation='none')
    plt.colorbar(fig2, shrink=0.5)

    # Plot the Resid between them
    f3 = fig.add_subplot(gs[0,2])
    plt.title('Resid')
    fig3 = f3.imshow( (image - target_image).array,interpolation='none')
    plt.colorbar(fig3, shrink=0.5)
    plt.show()



    # Initialize the geometry of the grid
    gs = gridspec.GridSpec(1,3)
    fig = plt.figure(figsize=(10,5))

    # Plot the target data img
    f1 = fig.add_subplot(gs[0,0])
    plt.title('Target')
    fig1 = f1.imshow(target_image.array,interpolation='none',origin='lower')
    plt.colorbar(fig1, shrink=0.5) 

    # Plot the best fit img
    f2 = fig.add_subplot(gs[0,1])
    plt.title('Best Fit')
    fig2 = f2.imshow(image.array,interpolation='none')
    plt.colorbar(fig2, shrink=0.5)

    # Plot the Resid between them
    f3 = fig.add_subplot(gs[0,2])
    plt.title('Resid')
    fig3 = f3.imshow( (image - target_image).array,interpolation='none')
    plt.colorbar(fig3, shrink=0.5)
    plt.show()
    '''

#    ipdb.set_trace()
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
    
