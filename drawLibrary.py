# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 14:50:00 2014

@author: luis
"""

## Library containing sets of functions for drawing galaxies and residuals
## for usage with fitters.

from __future__ import division
import galsim
import numpy as np
import gauss
import lmfit
import scipy.optimize

# Draw from a distribution and return a binned image object with one galaxy. 
def drawShoot_galaxy(flux, hlr, e1, e2, x0, y0, x_len, y_len, scale, func, seed):
    big_fft_params = galsim.GSParams(maximum_fft_size=100240)    
    if func is galsim.Gaussian:    
        gal = func(half_light_radius=hlr, flux=flux, gsparams=big_fft_params)
        gal = gal.shear(g1=e1, g2=e2)
        gal = gal.shift(x0,y0)
        image = galsim.ImageD(x_len, y_len, scale=scale)
        image = gal.drawImage(image=image, method='phot',rng=seed)
        return image

# Use the analytic definition of an image profile for one galaxy. 
def draw_galaxy_1(flux, hlr, e1, e2, x0, y0, x_len, y_len, scale, func):
    big_fft_params = galsim.GSParams(maximum_fft_size=100240)
    if func is galsim.Gaussian:        
        gal = func(half_light_radius=hlr, flux=flux, gsparams=big_fft_params)
        gal = gal.shear(g1=e1, g2=e2)
        gal = gal.shift(x0,y0)
        image = galsim.ImageD(x_len, y_len, scale=scale)
        image = gal.drawImage(image=image)
        return image

# The difference of the data and the model for one galaxy
# that is to be reduced.
def resid_1(param, target_image, x_len, y_len, scale):
    flux = param['flux'].value
    hlr = param['hlr'].value
    e1 = param['e1'].value
    e2 = param['e2'].value
    x0 = param['x0'].value
    y0 = param['y0'].value
    image = draw_galaxy_1(flux,hlr,e1,e2,x0,y0,x_len,y_len,scale)
    return (image-target_image).array.ravel()

# Draw from two distributions and return a binned image object with two 
# galaxies. 
def drawShoot_galaxy_2(flux_1,hlr_1,e1_1,e2_1,x_center1,y_center1,
                       flux_2,hlr_2,e1_2,e2_2,x_center2,y_center2,
                       x_len,y_len,scale,func_1,func_2,seed):
    big_fft_params = galsim.GSParams(maximum_fft_size=100240)                       
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

# Use the analytic definition of an image profile for two galaxies. 
def draw_galaxy_2(flux_1,hlr_1,e1_1,e2_1,x_center1,y_center1,
                  flux_2,hlr_2,e1_2,e2_2,x_center2,y_center2,
                  x_len,y_len,scale,func_1,func_2):
    big_fft_params = galsim.GSParams(maximum_fft_size=100240)
    if func_1 is galsim.Gaussian:                       
        gal_1 = func_1(half_light_radius=hlr_1, flux=flux_1, gsparams=big_fft_params)
        gal_1 = gal_1.shear(g1=e1_1, g2=e2_1)
        gal_1 = gal_1.shift(x_center1,y_center1)
        image_1 = galsim.ImageD(x_len, y_len, scale=scale)
        image_1 = gal_1.drawImage(image=image_1)
    if func_2 is galsim.Gaussian:
        gal_2 = func_2(half_light_radius=hlr_2, flux=flux_2, gsparams=big_fft_params)
        gal_2 = gal_2.shear(g1=e1_2, g2=e2_2)
        gal_2 = gal_2.shift(x_center2,y_center2)
        image_2 = galsim.ImageD(x_len, y_len, scale=scale)
        image_2 = gal_2.drawImage(image=image_2)
    image = image_1 + image_2
    
    return image

# The difference of the data and the model for one galaxy
# that is to be reduced.
def resid_2(param, target_image, x_len, y_len, scale, func_1, func_2):
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

# Function definition to return the original data array, best-fit array,
# residual, and correlation matrix with differences and error on e1 and e2.
def run_2_galaxy_full_params(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                             flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                             size_1,size_2,pixel_scale,func_gauss_1,func_gauss_2,dev_1,dev_2):

    im_1 = drawShoot_galaxy(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                        size_1,size_2,pixel_scale,func_gauss_1,dev_1)
    im_2 = drawShoot_galaxy(flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                                        size_1,size_2,pixel_scale,func_gauss_2,dev_2)
    
    im = im_1 + im_2
    
    # Obtain the image bounds and domain information
    x_cen,y_cen,x,y,X,Y = return_domain(im)
    
    # -----------------------------------------------------------------------    
    # Estimate the parameters of the image.
    
    # Define some seed that's far from true values and insert into
    # lmfit object for galaxy one and two
    p0 = 1.0*np.array([flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                       flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b])
    parameters = lmfit.Parameters()
    parameters.add('flux_1', value=p0[0])
    parameters.add('hlr_1', value=p0[1], min=0.0)
    parameters.add('e1_1', value=p0[2], min=-1.0, max=1.0)
    parameters.add('e2_1', value=p0[3], min=-1.0, max=1.0)
    parameters.add('x_center1',value=p0[4])
    parameters.add('y_center1',value=p0[5])
    
    parameters.add('flux_2', value=p0[6])
    parameters.add('hlr_2', value=p0[7], min=0.0)
    parameters.add('e1_2', value=p0[8], min=-1.0, max=1.0)
    parameters.add('e2_2', value=p0[9], min=-1.0, max=1.0)
    parameters.add('x_center2',value=p0[10])
    parameters.add('y_center2',value=p0[11])
    
    
    # Extract params that minimize the difference of the data from the model.
    result = lmfit.minimize(resid_2, parameters, args=(im,size_1,size_2,pixel_scale,func_gauss_1,func_gauss_2))
    best_fit = draw_galaxy_2(result.params['flux_1'].value,
                                       result.params['hlr_1'].value,
                                       result.params['e1_1'].value,
                                       result.params['e2_1'].value,
                                       result.params['x_center1'].value,
                                       result.params['y_center1'].value,
                                       result.params['flux_2'].value,
                                       result.params['hlr_2'].value,
                                       result.params['e1_2'].value,
                                       result.params['e2_2'].value,
                                       result.params['x_center2'].value,
                                       result.params['y_center2'].value,
                                       size_1,size_2,pixel_scale,func_gauss_1,func_gauss_2)
                                                                      
    return im,best_fit,result
        

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
    
    
## New Model for Parameterization--------------------------------------------

# Create a two-system model. (Will change parameters with time)
def model_objects_a_b_gaussians(sum_flux,dif_flux,x_sum,x_dif,y_sum,y_dif,
                                hlr_a,e1_a,e2_a,
                                hlr_b,e1_b,e2_b,
                                x_len,y_len,scale,func_a,func_b):
                                    
                                    
    big_fft_params = galsim.GSParams(maximum_fft_size=100240)    
    
    # Calculate the individual flux's from the sum and the difference.
    flux_a = (sum_flux + dif_flux)/2.0; flux_b = (sum_flux - dif_flux)/2.0
    # Calculate the coordinates of the centroids
    x0_a = (x_sum + x_dif)*(flux_a + flux_b)/(2*flux_a)
    x0_b = (x_sum - x_dif)*(flux_a + flux_b)/(2*flux_b)
    y0_a = (y_sum + y_dif)*(flux_a + flux_b)/(2*flux_a)
    y0_b = (y_sum - y_dif)*(flux_a + flux_b)/(2*flux_b)
    
    if func_a is galsim.Gaussian:                       
        gal_1 = func_a(half_light_radius=hlr_a, flux=flux_a, gsparams=big_fft_params)
        gal_1 = gal_1.shear(e1=e1_a, e2=e2_a)
        gal_1 = gal_1.shift(x0_a,y0_a)
        image_1 = galsim.ImageD(x_len, y_len, scale=scale)
        image_1 = gal_1.drawImage(image=image_1)
    if func_b is galsim.Gaussian:
        gal_2 = func_b(half_light_radius=hlr_b, flux=flux_b, gsparams=big_fft_params)
        gal_2 = gal_2.shear(e1=e1_b, e2=e2_b)
        gal_2 = gal_2.shift(x0_b,y0_b)
        image_2 = galsim.ImageD(x_len, y_len, scale=scale)
        image_2 = gal_2.drawImage(image=image_2)
    image = image_1 + image_2
    
    return image
    
def resid_objects_a_b_gaussians(param,data_image,x_len,y_len,scale,func_a,func_b):
    sum_flux = param['sum_flux'].value
    dif_flux = param['dif_flux'].value    
    
    x_sum = param['x_sum'].value
    x_dif = param['x_dif'].value    
    y_sum = param['y_sum'].value
    y_dif = param['y_dif'].value    
    
    hlr_a = param['hlr_a'].value
    e1_a = param['e1_a'].value
    e2_a = param['e2_a'].value
    
    hlr_b = param['hlr_b'].value
    e1_b = param['e1_b'].value
    e2_b = param['e2_b'].value

    image = model_objects_a_b_gaussians(sum_flux,dif_flux,x_sum,x_dif,y_sum,y_dif,
                                        hlr_a,e1_a,e2_a,
                                        hlr_b,e1_b,e2_b,
                                        x_len,y_len,scale,func_a,func_b)    
        
    return (image - data_image).array.ravel()

    
def run_2_galaxy(sum_flux,dif_flux,x_sum,x_dif,y_sum,y_dif,
                 hlr_a,e1_a,e2_a,
                 hlr_b,e1_b,e2_b,
                 x_len,y_len,pixel_scale,func_gauss_1,func_gauss_2,dev_1,dev_2,dev_3):

    # Define the variables needed to produce an image    
    flux_a = (sum_flux + dif_flux)/2.0; flux_b = (sum_flux - dif_flux)/2.0
    
    # Calculate the coordinates of the centroids
    x0_a = (x_sum + x_dif)*(flux_a + flux_b)/(2*flux_a)
    x0_b = (x_sum - x_dif)*(flux_a + flux_b)/(2*flux_b)
    y0_a = (y_sum + y_dif)*(flux_a + flux_b)/(2*flux_a)
    y0_b = (y_sum - y_dif)*(flux_a + flux_b)/(2*flux_b)
        
    # Create the image
    im_1 = drawShoot_galaxy(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,
                            x_len,y_len,pixel_scale,func_gauss_1,dev_1)
    im_2 = drawShoot_galaxy(flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,
                            x_len,y_len,pixel_scale,func_gauss_2,dev_2)
    

    im = im_1 + im_2
    #im.addNoise(galsim.PoissonNoise(rng=dev_3,sky_level=sky_level))
    
     # Obtain the image bounds and domain information
    x_cen,y_cen,x,y,X,Y = return_domain(im)
    
    # -----------------------------------------------------------------------    
    # Estimate the parameters of the image.
    
    # Define some seed that's far from true values and insert into
    # lmfit object for galaxy a and b
    p0 = 1.0*np.array([sum_flux,dif_flux,x_sum,x_dif,y_sum,y_dif,
                       hlr_a,e1_a,e2_a,
                       hlr_b,e1_b,e2_b])
    
    # Add the parameters
    parameters = lmfit.Parameters()
    parameters.add('sum_flux', value=p0[0])
    parameters.add('dif_flux', value=p0[1])
    
    parameters.add('x_sum',value=p0[2])
    parameters.add('x_dif',value=p0[3])    
    
    parameters.add('y_sum',value=p0[4])
    parameters.add('y_dif',value=p0[5])    
    
    parameters.add('hlr_a', value=p0[6], min=0.0)
    parameters.add('e1_a', value=p0[7], min=-1.0, max=1.0)
    parameters.add('e2_a', value=p0[8], min=-1.0, max=1.0)
    
    parameters.add('hlr_b', value=p0[9], min=0.0)
    parameters.add('e1_b', value=p0[10], min=-1.0, max=1.0)
    parameters.add('e2_b', value=p0[11], min=-1.0, max=1.0)
    
    # Extract params that minimize the difference of the data from the model.
    result = lmfit.minimize(resid_objects_a_b_gaussians, parameters, args=(im,x_len,y_len,pixel_scale,func_gauss_1,func_gauss_2))
    best_fit = model_objects_a_b_gaussians(result.params['sum_flux'].value,
                                           result.params['dif_flux'].value,
                                           result.params['x_sum'].value,
                                           result.params['x_dif'].value,
                                           
                                           result.params['y_sum'].value,
                                           result.params['y_dif'].value,                                           
                                           
                                           result.params['hlr_a'].value,                                            
                                           result.params['e1_a'].value,
                                           result.params['e2_a'].value,
                                           
                                           result.params['hlr_b'].value,
                                           result.params['e1_b'].value,
                                           result.params['e2_b'].value,
                                        
                                           x_len,y_len,pixel_scale,func_gauss_1,func_gauss_2)
                                           
    # Return the original image, the best fit, and the result of the fit.                                       
    return im,best_fit,result    

# Return the second central moment matrix from the Half-light radius,
# e1, and e2.     
def Q_Gauss(hlr,e1,e2):
    """ 
    Using Mathematica-generated fitting function, compute the circularized second moment
    square radius sqrt(r^2) = sqrt(Ixx + Iyy); where Ixx = Iyy and Ixy = 0.0
    """
    n=1/2
    r = hlr * (0.985444 + n * (0.391016 + n * (0.0739602
        + n * (0.00698719 + n * (0.00212432
        + n * (-0.000154052 + n * 0.0000219632))))))
    
    q = (1-np.sqrt(e1**2 + e2**2))/(1+np.sqrt(e1**2 + e2**2))
                                
    unrotI = 0.5 * r**2 * np.matrix([[1./q, 0.0],
                                      [0.0, q]], dtype=float)

    if e1 == 0 and e2 == 0:    
        phi = 0
    elif e1 == 0 and e2 != 0:
        phi = np.pi/2
    else:    
        phi = np.arctan2(e2,e1)                                      
                                      
    cph = np.cos(phi)
    sph = np.sin(phi)
    R = np.matrix([[cph, -sph],
                   [sph, cph]], dtype=float)
    return R * unrotI * R.T

# Compute the 2nd central moments of the mixture
def Q_two_mixture(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,
                  flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b):

    Q_a = Q_Gauss(hlr_a,e1_a,e2_a)
    Q_b = Q_Gauss(hlr_b,e1_b,e2_b)
    
    flux_tot = flux_a + flux_b
    
    mu_a = np.array([x0_a,y0_a])
    mu_b = np.array([x0_b,y0_b])
    weighted_mu = (flux_a*mu_a + flux_b*mu_b)/flux_tot
    
    return (1/flux_tot)*(flux_a*Q_a + flux_b*Q_b + 
            flux_a*np.outer(mu_a-weighted_mu,mu_a-weighted_mu) + 
            flux_b*np.outer(mu_b-weighted_mu,mu_b-weighted_mu))
            
# Compute the 2nd central moments of the mixture using non-linear solver
def Q_Gauss_solver(hlr,e1,e2):
    def F(x):
        n=1/2
        r = hlr * (0.985444 + n * (0.391016 + n * (0.0739602
        + n * (0.00698719 + n * (0.00212432
        + n * (-0.000154052 + n * 0.0000219632))))))
        return np.array([(x[0]-x[2])/(x[0]+x[2]+2*(x[0]*x[2]-x[1]**2)**(1/2)),2*x[1]/(x[0]+x[2]+2*(x[0]*x[2]-x[1]**2)**(1/2)),x[0]+x[2]]) - np.array([e1,e2,r**2])
               
    soln = scipy.optimize.broyden1(F,[.2,.1,.5],f_tol=1e-14)
    return np.mat([soln[0],soln[2]],
                  [soln[2],soln[1]])   
                  
# Extract the physical parameters from the result
def extract_physical_params(result):
    sum_flux = result.params['sum_flux'].value
    dif_flux = result.params['dif_flux'].value
    x_sum = result.params['x_sum'].value
    x_dif = result.params['x_dif'].value
       
    y_sum = result.params['y_sum'].value
    y_dif = result.params['y_dif'].value                                           
       
    hlr_a = result.params['hlr_a'].value                                            
    e1_a = result.params['e1_a'].value
    e2_a = result.params['e2_a'].value
       
    hlr_b = result.params['hlr_b'].value
    e1_b = result.params['e1_b'].value
    e2_b = result.params['e2_b'].value
    
    flux_a = (sum_flux + dif_flux)/2.0; flux_b = (sum_flux - dif_flux)/2.0
    
    # Calculate the coordinates of the centroids
    x0_a = (x_sum + x_dif)*(flux_a + flux_b)/(2*flux_a)
    x0_b = (x_sum - x_dif)*(flux_a + flux_b)/(2*flux_b)
    y0_a = (y_sum + y_dif)*(flux_a + flux_b)/(2*flux_a)
    y0_b = (y_sum - y_dif)*(flux_a + flux_b)/(2*flux_b)
    
    return flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b
    
    