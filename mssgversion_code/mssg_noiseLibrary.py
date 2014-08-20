# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 13:57:03 2014

@author: luis
"""

# mssg version 
# changes:
#  --- added fit_2_gals function


## Library containing sets of functions for drawing galaxies and residuals
## with a modified code structure for optional arguments, psf, and noise.

from __future__ import division
import galsim
import numpy as np
import gauss
import lmfit
import scipy.optimize
import drawLibrary
import matplotlib.pyplot as plt
import matplotlib.gridspec as grid 

# Create a galaxy with a sersic profile and optional psf to the image. 
def create_galaxy(flux, hlr, e1, e2, x0, y0, galtype_gal=galsim.Sersic, sersic_index=0.5,
                  psf_flag=False, psf_type=galsim.Moffat, beta=5, size_psf=1, flux_psf=1,
                  x_len=100, y_len=100, scale=0.2, method='fft',seed=None):
    big_fft_params = galsim.GSParams(maximum_fft_size=1002400)

#    print "\nPostage Stamp is", x_len, "by", y_len, "with\na scale of", scale,"\"/Pixel"    
        
    if galtype_gal is galsim.Sersic:
        assert sersic_index != 0
        if sersic_index == 0.5:
            print "\nThe object drawn is a gaussian with n = 0.5" 
        gal = galtype_gal(n=sersic_index, half_light_radius=hlr, flux=flux, gsparams=big_fft_params)
        gal = gal.shear(g1=e1, g2=e2)
        gal = gal.shift(x0,y0)
        image = galsim.ImageD(x_len, y_len, scale=scale)
        if psf_flag == True:
            psf_gal = convolve_with_psf(gal, beta=beta, size_psf=size_psf, psf_type=psf_type, flux_psf=flux_psf)
            if method == 'fft':
                image = psf_gal.drawImage(image=image, method=method)
            else:
                image = psf_gal.drawImage(image=image, method=method,rng=seed)
            return image
        else:
            if method == 'fft':
                image = gal.drawImage(image=image, method=method)
            else:
                image = gal.drawImage(image=image, method=method,rng=seed)
            return image    
        
    else:
        raise ValueError("Not using a sersic profile for the object.")
        
# Draw objects onto an image or return the objects themselves.
def draw_simple(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,
                flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,n_b,
                psf_flag,beta,fwhm_psf,
                x_len,y_len,pixel_scale,galtype_a,galtype_b,seed_a,seed_b,seed_p,
                add_noise_flag,sky_level):

    image_a = create_galaxy(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,galtype_gal=galtype_a,sersic_index=n_a,
                            x_len=x_len,y_len=y_len,scale=pixel_scale,
                            psf_flag=psf_flag, beta=beta, size_psf=fwhm_psf)
                                
    image_b = create_galaxy(flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,galtype_gal=galtype_b,sersic_index=n_b,
                            x_len=x_len,y_len=y_len,scale=pixel_scale,
                            psf_flag=psf_flag, beta=beta, size_psf=fwhm_psf)
                            
    image = image_a + image_b
    if add_noise_flag == True:
        image = add_noise(image,seed=seed_p,sky_level=sky_level)        
    return image
    
# Add poisson noise to an image with a sky level argument.        
def add_noise(image, noise_type=galsim.PoissonNoise, seed=None, sky_level=0):
    if sky_level == 0: 
        print "Adding Poisson noise without sky level."     
    else: 
        print "Adding Poisson noise with sky level = ", sky_level
    if noise_type is galsim.PoissonNoise:    
        image.addNoise(noise_type(sky_level=sky_level,rng=seed))
        return image
    else:
        raise ValueError("Not using poisson noise in your image.")
        
def calc_SNR(im, texp, sbar, weight):
    # Work with the image in count per second
    Di = im.array/texp
    # Calculate the threshold and mask per second
    threshold_per_s = weight*np.sqrt(sbar/texp)
    mask_per_s = Di > threshold_per_s
    # Calculate SNR using David's formulation
    nu_s = np.sqrt(texp/sbar)*np.sqrt((mask_per_s*Di*Di).sum())
    
    # Now work with the original galsim image (Di*texp=im.array)
    threshold = weight*np.sqrt(texp*sbar)
    mask = im.array > threshold
    # Now calculate the SNR using the original galsim image
    nu = np.sqrt(1/(sbar*texp))*np.sqrt((mask*im.array**2).sum())

    return nu_s, mask_per_s, nu, mask

# Convolve an object with a PSF.
def convolve_with_psf(gal, beta, size_psf, psf_type=galsim.Moffat, flux_psf=1):
    big_fft_params = galsim.GSParams(maximum_fft_size=1002400)
    print "Using a psf with beta =", beta,"and size = ", size_psf," \"" 
    psf = psf_type(beta=beta, fwhm=size_psf, flux=flux_psf, gsparams=big_fft_params)
    psf_gal = galsim.Convolve([gal,psf])
    return psf_gal
    
    
def residual_func_simple(param, data_image, sky_level, x_len, y_len, pixel_scale, 
                         galtype_a,n_a,galtype_b,n_b,galtype,n):
    # If we are looking at one object only:
    if len(param) == 6:
        assert galtype != None
        flux = param['flux'].value
        hlr = param['hlr'].value
        e1 = param['e1'].value
        e2 = param['e2'].value
        x0 = param['x0'].value
        y0 = param['y0'].value
        image = create_galaxy(flux,hlr,e1,e2,x0,y0,galtype_gal=galtype,sersic_index=n,
                              x_len=x_len,y_len=y_len,scale=pixel_scale)
        # Create error array
        error = np.sqrt(data_image.array.ravel())
        # Set the errors equal to 1 where 0 is.
        error[error==0] = 1
        return (image-data_image).array.ravel()
        
    # Else we are looking at two:    
    else:
        assert galtype_a != None
        assert galtype_b != None
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
        image_a = create_galaxy(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,galtype_gal=galtype_a,sersic_index=n_a,
                                x_len=x_len,y_len=y_len,scale=pixel_scale)
                                
        image_b = create_galaxy(flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,galtype_gal=galtype_b,sersic_index=n_b,
                                x_len=x_len,y_len=y_len,scale=pixel_scale)
                            
        image = image_a + image_b
        
        if sky_level > 10:        
            return (image-data_image).array.ravel()/np.sqrt(sky_level + data_image.array).ravel()
        else:
            return (image-data_image).array.ravel()

        
# Function definition to return the original data array, best-fit array,
# residual, and correlation matrix with differences and error on e1 and e2.
def run_2_galaxy_full_params_simple(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,
                                    flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,n_b,
                                    psf_flag,beta,fwhm_psf,
                                    x_len,y_len,pixel_scale,galtype_a,galtype_b,seed_a,seed_b,seed_p,
                                    add_noise_flag,sky_level):

    image_a = create_galaxy(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,galtype_gal=galtype_a,sersic_index=n_a,
                            x_len=x_len,y_len=y_len,scale=pixel_scale,
                            psf_flag=psf_flag, beta=beta, size_psf=fwhm_psf)
                                
    image_b = create_galaxy(flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,galtype_gal=galtype_b,sersic_index=n_b,
                            x_len=x_len,y_len=y_len,scale=pixel_scale,
                            psf_flag=psf_flag, beta=beta, size_psf=fwhm_psf)
                            
    image = image_a + image_b
    if add_noise_flag == True:
        image = add_noise(image,seed=seed_p,sky_level=sky_level)    
    # Obtain the image bounds and domain information
    x_cen,y_cen,x,y,X,Y = drawLibrary.return_domain(image)
    
    # -----------------------------------------------------------------------    
    # Estimate the parameters of the image.
    
    # Define some seed that's far from true values and insert into
    # lmfit object for galaxy one and two
    p0 = 1.0*np.array([flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,
                       flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b])
    parameters = lmfit.Parameters()
    parameters.add('flux_a', value=p0[0])
    parameters.add('hlr_a', value=p0[1], min=0.0)
    parameters.add('e1_a', value=p0[2], min=-1.0, max=1.0)
    parameters.add('e2_a', value=p0[3], min=-1.0, max=1.0)
    parameters.add('x0_a',value=p0[4])
    parameters.add('y0_a',value=p0[5])
    
    parameters.add('flux_b', value=p0[6])
    parameters.add('hlr_b', value=p0[7], min=0.0)
    parameters.add('e1_b', value=p0[8], min=-1.0, max=1.0)
    parameters.add('e2_b', value=p0[9], min=-1.0, max=1.0)
    parameters.add('x0_b',value=p0[10])
    parameters.add('y0_b',value=p0[11])
    
    
    # Extract params that minimize the difference of the data from the model.
    result = lmfit.minimize(residual_func_simple, parameters, args=(image, sky_level, x_len, y_len, pixel_scale, galtype_a, n_a, galtype_b, n_b, None,0))
    best_fit_a = create_galaxy(result.params['flux_a'].value,
                               result.params['hlr_a'].value,
                               result.params['e1_a'].value,
                               result.params['e2_a'].value,
                               result.params['x0_a'].value,
                               result.params['y0_a'].value,
                               x_len=x_len,y_len=y_len,scale=pixel_scale,galtype_gal=galtype_a,sersic_index=n_a)
                               
    best_fit_b = create_galaxy(result.params['flux_b'].value,
                               result.params['hlr_b'].value,
                               result.params['e1_b'].value,
                               result.params['e2_b'].value,
                               result.params['x0_b'].value,
                               result.params['y0_b'].value,
                               x_len=x_len,y_len=y_len,scale=pixel_scale,galtype_gal=galtype_b,sersic_index=n_b)
                                       
                                                                      
    return image, best_fit_a+best_fit_b, result


def fit_2_galaxies(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,
                                    flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,n_b,
                                    psf_flag,beta,fwhm_psf,
                                    x_len,y_len,pixel_scale,galtype_a,galtype_b,seed_a,seed_b,seed_p,
                                    add_noise_flag,sky_level, galimg ):
    '''
    image_a = create_galaxy(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,galtype_gal=galtype_a,sersic_index=n_a,
                            x_len=x_len,y_len=y_len,scale=pixel_scale,
                            psf_flag=psf_flag, beta=beta, size_psf=fwhm_psf)
                                
    image_b = create_galaxy(flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,galtype_gal=galtype_b,sersic_index=n_b,
                            x_len=x_len,y_len=y_len,scale=pixel_scale,
                            psf_flag=psf_flag, beta=beta, size_psf=fwhm_psf)
                            
    image = image_a + image_b

    if add_noise_flag == True:
        image = add_noise(image,seed=seed_p,sky_level=sky_level)    
        '''


    image = galimg

    # Obtain the image bounds and domain information
    x_cen,y_cen,x,y,X,Y = drawLibrary.return_domain(image)
    
    # -----------------------------------------------------------------------    
    # Estimate the parameters of the image.
    
    # Define some seed that's far from true values and insert into
    # lmfit object for galaxy one and two
    p0 = 1.0*np.array([flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,
                       flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b])
    parameters = lmfit.Parameters()
    parameters.add('flux_a', value=p0[0])
    parameters.add('hlr_a', value=p0[1], min=0.0)
    parameters.add('e1_a', value=p0[2], min=-1.0, max=1.0)
    parameters.add('e2_a', value=p0[3], min=-1.0, max=1.0)
    parameters.add('x0_a',value=p0[4])
    parameters.add('y0_a',value=p0[5])
    
    parameters.add('flux_b', value=p0[6])
    parameters.add('hlr_b', value=p0[7], min=0.0)
    parameters.add('e1_b', value=p0[8], min=-1.0, max=1.0)
    parameters.add('e2_b', value=p0[9], min=-1.0, max=1.0)
    parameters.add('x0_b',value=p0[10])
    parameters.add('y0_b',value=p0[11])
    
    
    # Extract params that minimize the difference of the data from the model.
    result = lmfit.minimize(residual_func_simple, parameters, args=(image, sky_level, x_len, y_len, pixel_scale, galtype_a, n_a, galtype_b, n_b, None,0))
    best_fit_a = create_galaxy(result.params['flux_a'].value,
                               result.params['hlr_a'].value,
                               result.params['e1_a'].value,
                               result.params['e2_a'].value,
                               result.params['x0_a'].value,
                               result.params['y0_a'].value,
                               x_len=x_len,y_len=y_len,scale=pixel_scale,galtype_gal=galtype_a,sersic_index=n_a)
                               
    best_fit_b = create_galaxy(result.params['flux_b'].value,
                               result.params['hlr_b'].value,
                               result.params['e1_b'].value,
                               result.params['e2_b'].value,
                               result.params['x0_b'].value,
                               result.params['y0_b'].value,
                               x_len=x_len,y_len=y_len,scale=pixel_scale,galtype_gal=galtype_b,sersic_index=n_b)
                                       
                                                                      
    return image, best_fit_a+best_fit_b, result

def residual_func_complex(param, data_image, sky_level, x_len, y_len, pixel_scale, 
                          galtype_a_bulge, galtype_a_disk, n_a_bulge, n_a_disk,
                          galtype_b_bulge, galtype_b_disk, n_b_bulge, n_b_disk,
                          galtype_bulge, galtype_disk, n_bulge, n_disk):
    # If we are looking at one object only:
    if len(param) == 6:
        assert galtype_bulge != None
        assert galtype_disk != None
        flux_tot = param['flux_tot'].value
        bulge_frac = param['bulge_frac'].value        
        hlr_bulge = param['hlr_bulge'].value
        hlr_disk = param['hlr_disk'].value
        e1 = param['e1'].value
        e2 = param['e2'].value
        x0 = param['x0'].value
        y0 = param['y0'].value
        flux_bulge = flux_tot*bulge_frac
        flux_disk = flux_tot*(1-bulge_frac)
        image_bulge = create_galaxy(flux_bulge,hlr_bulge,e1,e2,x0,y0,galtype_gal=galtype_bulge,sersic_index=n_bulge,
                                    x_len=x_len,y_len=y_len,scale=pixel_scale)
        image_disk = create_galaxy(flux_disk,hlr_disk,e1,e2,x0,y0,galtype_gal=galtype_disk,sersic_index=n_disk,
                                    x_len=x_len,y_len=y_len,scale=pixel_scale) 
        image = image_bulge + image_disk                            
        # Create error array
        error = np.sqrt(data_image.array.ravel())
        # Set the errors equal to 1 where 0 is.
        error[error==0] = 1
        return (image-data_image).array.ravel()
        
    # Else we are looking at two with one component each:    
    else:
        assert galtype_a_bulge != None
        assert galtype_a_disk != None
        assert galtype_b_bulge != None
        assert galtype_b_disk != None
        flux_a_tot = param['flux_a_tot'].value
        bulge_a_frac = param['bulge_a_frac'].value        
        hlr_a_bulge = param['hlr_a_bulge'].value
        hlr_a_disk = param['hlr_a_disk'].value        
        e1_a = param['e1_a'].value
        e2_a = param['e2_a'].value
        x0_a = param['x0_a'].value
        y0_a = param['y0_a'].value
        # Extract meaningful params
        flux_a_bulge = flux_a_tot*bulge_a_frac
        flux_a_disk = flux_a_tot*(1-bulge_a_frac)

        flux_b_tot = param['flux_b_tot'].value
        bulge_b_frac = param['bulge_b_frac'].value        
        hlr_b_bulge = param['hlr_b_bulge'].value
        hlr_b_disk = param['hlr_b_disk'].value            
        e1_b = param['e1_b'].value
        e2_b = param['e2_b'].value
        x0_b = param['x0_b'].value
        y0_b = param['y0_b'].value
        # Extract meaningful params
        flux_b_bulge = flux_b_tot*bulge_b_frac
        flux_b_disk = flux_b_tot*(1-bulge_b_frac)

        
        image_a_bulge = create_galaxy(flux_a_bulge,hlr_a_bulge,e1_a,e2_a,x0_a,y0_a,galtype_gal=galtype_a_bulge,sersic_index=n_a_bulge,
                                      x_len=x_len,y_len=y_len,scale=pixel_scale)
        image_a_disk = create_galaxy(flux_a_disk,hlr_a_disk,e1_a,e2_a,x0_a,y0_a,galtype_gal=galtype_a_disk,sersic_index=n_a_disk,
                                     x_len=x_len,y_len=y_len,scale=pixel_scale)
    

                                
        image_b_bulge = create_galaxy(flux_b_bulge,hlr_b_bulge,e1_b,e2_b,x0_b,y0_b,galtype_gal=galtype_b_bulge,sersic_index=n_b_bulge,
                                      x_len=x_len,y_len=y_len,scale=pixel_scale)

        image_b_disk = create_galaxy(flux_b_disk,hlr_b_disk,e1_b,e2_b,x0_b,y0_b,galtype_gal=galtype_b_bulge,sersic_index=n_b_disk,
                                     x_len=x_len,y_len=y_len,scale=pixel_scale)

                            
        image = image_a_bulge + image_a_disk + image_b_bulge + image_b_disk
        
        if sky_level > 10:        
            return (image-data_image).array.ravel()
        else:
            return (image-data_image).array.ravel()
        
# Function definition to return the original data array, best-fit array,
# residual, and correlation matrix with differences and error on e1 and e2.
def run_2_galaxy_full_params_complex(flux_a_tot,bulge_a_frac,hlr_a_bulge,n_a_bulge,hlr_a_disk,n_a_disk,e1_a,e2_a,x0_a,y0_a,
                                     flux_b_tot,bulge_b_frac,hlr_b_bulge,n_b_bulge,hlr_b_disk,n_b_disk,e1_b,e2_b,x0_b,y0_b,
                                     psf_flag,beta,fwhm_psf,
                                     x_len,y_len,pixel_scale,
                                     galtype_a,galtype_b,
                                     seed_a,seed_b,seed_p,
                                     add_noise_flag,sky_level):

    flux_a_bulge = flux_a_tot*bulge_a_frac
    flux_a_disk = flux_a_tot*(1-bulge_a_frac)
    
    image_a_bulge = create_galaxy(flux_a_bulge,hlr_a_bulge,e1_a,e2_a,x0_a,y0_a,galtype_gal=galtype_a,sersic_index=n_a_bulge,
                                  x_len=x_len,y_len=y_len,scale=pixel_scale,
                                  psf_flag=psf_flag, beta=beta, size_psf=fwhm_psf)
                                  
    image_a_disk = create_galaxy(flux_a_disk,hlr_a_disk,e1_a,e2_a,x0_a,y0_a,galtype_gal=galtype_a,sersic_index=n_a_disk,
                                  x_len=x_len,y_len=y_len,scale=pixel_scale,
                                  psf_flag=psf_flag, beta=beta, size_psf=fwhm_psf)
    
    flux_b_bulge = flux_b_tot*bulge_b_frac
    flux_b_disk = flux_b_tot*(1-bulge_b_frac)
                                
    image_b_bulge = create_galaxy(flux_b_bulge,hlr_b_bulge,e1_b,e2_b,x0_b,y0_b,galtype_gal=galtype_b,sersic_index=n_b_bulge,
                                  x_len=x_len,y_len=y_len,scale=pixel_scale,
                                  psf_flag=psf_flag, beta=beta, size_psf=fwhm_psf)

    image_b_disk = create_galaxy(flux_b_disk,hlr_b_disk,e1_b,e2_b,x0_b,y0_b,galtype_gal=galtype_b,sersic_index=n_b_disk,
                                  x_len=x_len,y_len=y_len,scale=pixel_scale,
                                  psf_flag=psf_flag, beta=beta, size_psf=fwhm_psf)

                            
    image = image_a_bulge + image_a_disk + image_b_bulge + image_b_disk
    if add_noise_flag == True:
        image = add_noise(image,seed=seed_p,sky_level=sky_level)    
    # Obtain the image bounds and domain information
    x_cen,y_cen,x,y,X,Y = drawLibrary.return_domain(image)
    
    # -----------------------------------------------------------------------    
    # Estimate the parameters of the image.
    
    # Define some seed that's far from true values and insert into
    # lmfit object for galaxy one and two
    p0 = 1.0*np.array([flux_a_tot,bulge_a_frac,hlr_a_bulge,hlr_a_disk,e1_a,e2_a,x0_a,y0_a,
                       flux_b_tot,bulge_b_frac,hlr_b_bulge,hlr_b_disk,e1_b,e2_b,x0_b,y0_b])
    parameters = lmfit.Parameters()
    parameters.add('flux_a_tot', value=p0[0])
    parameters.add('bulge_a_frac',value=p0[1])
    parameters.add('hlr_a_bulge', value=p0[2])
    parameters.add('hlr_a_disk',value=p0[3])
    parameters.add('e1_a', value=p0[4], min=-1.0, max=1.0)
    parameters.add('e2_a', value=p0[5], min=-1.0, max=1.0)
    parameters.add('x0_a',value=p0[6])
    parameters.add('y0_a',value=p0[7])
    
    parameters.add('flux_b_tot', value=p0[8])
    parameters.add('bulge_b_frac',value=p0[9])
    parameters.add('hlr_b_bulge', value=p0[10])
    parameters.add('hlr_b_disk',value=p0[11])
    parameters.add('e1_b', value=p0[12], min=-1.0, max=1.0)
    parameters.add('e2_b', value=p0[13], min=-1.0, max=1.0)
    parameters.add('x0_b',value=p0[14])
    parameters.add('y0_b',value=p0[15])
    
    
    # Extract params that minimize the difference of the data from the model.
    result = lmfit.minimize(residual_func_complex, parameters, 
                            args=(image, sky_level, x_len, y_len, pixel_scale, 
                                  galtype_a, galtype_a, n_a_bulge, n_a_disk,
                                  galtype_b, galtype_b, n_b_bulge, n_b_disk,
                                  None,None,0,0)) # Case for one object
                                  
    # Extract the  params to draw the best fit objects.     
    flux_a_tot_est = result.params['flux_a_tot'].value
    bulge_a_frac_est = result.params['bulge_a_frac'].value
    hlr_a_bulge_est = result.params['hlr_a_bulge'].value
    hlr_a_disk_est = result.params['hlr_a_disk'].value
    e1_a_est = result.params['e1_a'].value
    e2_a_est = result.params['e2_a'].value
    x0_a_est = result.params['x0_a'].value
    y0_a_est = result.params['y0_a'].value
    # Extract meaningful params for drawing.
    flux_a_bulge_est = flux_a_tot_est*bulge_a_frac_est
    flux_a_disk_est = flux_a_tot_est*(1-bulge_a_frac_est)

    flux_b_tot_est = result.params['flux_b_tot'].value
    bulge_b_frac_est = result.params['bulge_b_frac'].value
    hlr_b_bulge_est = result.params['hlr_b_bulge'].value
    hlr_b_disk_est = result.params['hlr_b_disk'].value
    e1_b_est = result.params['e1_b'].value
    e2_b_est = result.params['e2_b'].value
    x0_b_est = result.params['x0_b'].value
    y0_b_est = result.params['y0_b'].value
    # Extract menaingful params for drawing.
    flux_b_bulge_est = flux_b_tot_est*bulge_b_frac_est
    flux_b_disk_est = flux_b_tot_est*(1-bulge_b_frac_est)                              
                                  
    best_fit_a_bulge = create_galaxy(flux_a_bulge_est,hlr_a_bulge_est,e1_a_est,e2_a_est,x0_a_est,y0_a_est,
                                     x_len=x_len,y_len=y_len,scale=pixel_scale,galtype_gal=galtype_a,sersic_index=n_a_bulge)
                                     
    best_fit_a_disk = create_galaxy(flux_a_disk_est,hlr_a_disk_est,e1_a_est,e2_a_est,x0_a_est,y0_a_est,
                                     x_len=x_len,y_len=y_len,scale=pixel_scale,galtype_gal=galtype_a,sersic_index=n_a_disk)                                 

    best_fit_b_bulge = create_galaxy(flux_b_bulge_est,hlr_b_bulge_est,e1_b_est,e2_b_est,x0_b_est,y0_b_est,
                                     x_len=x_len,y_len=y_len,scale=pixel_scale,galtype_gal=galtype_a,sersic_index=n_b_bulge)
                                     
    best_fit_b_disk = create_galaxy(flux_b_disk_est,hlr_b_disk_est,e1_b_est,e2_b_est,x0_b_est,y0_b_est,
                                     x_len=x_len,y_len=y_len,scale=pixel_scale,galtype_gal=galtype_a,sersic_index=n_b_disk)                                                                        
                                                                      
    return image, best_fit_a_bulge + best_fit_a_disk + best_fit_b_bulge + best_fit_b_disk, result

# TODO implement more than one datum
def plot(domain, data, figsize, pos, row, col, colors, markers, legend=None, fontsize=12, legend_fontsize=12, t_fontsize=15, 
         text_x=0, text_y=0, text_s='', text_fs=12,
         yerr=None, suptitle='', x_label='', x_lim=None, y_label=None, y_lim=None, scatter_flag=True, plot_flag=True):
    figure = plt.figure(figsize=figsize)
    if text_s != '':
        plt.text(text_x,text_y,text_s,fontsize=text_fs)         
    gs = grid.GridSpec(row,col)
    for i, key in enumerate(data):
        figure.add_subplot(gs[pos[i]])
        plt.xlabel(x_label,fontsize=fontsize)
        if y_label != None:
            plt.ylabel(y_label[i],fontsize=fontsize)
        plt.suptitle(suptitle,fontsize=t_fontsize)
        if y_lim != None:
            plt.ylim(y_lim[i])
        if x_lim != None:
            plt.xlim(x_lim)
        if scatter_flag == True:
            p = plt.scatter(domain,data[key],c=colors[i],marker=markers[i])
            if legend != None:
                plt.legend([p],[legend[i]],loc=1,prop={'size':legend_fontsize})
        if plot_flag == True and scatter_flag == True:
            plt.plot(domain,data[key],c=colors[i])
        else:
            if legend != None:
                plt.plot(domain,data[key],c=colors[i],label=legend[i])
            else:
                plt.plot(domain,data[key],c=colors[i])
        if yerr != None:
            plt.errorbar(domain,data[key],c=colors[i])
    return figure        
                        
