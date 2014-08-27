# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 11:06:01 2014

@author: luis
"""

from __future__ import division
from pandas import Series, DataFrame
import pandas as pd
import numpy as np
import drawLibrary
import noiseLibrary
import galsim
import lmfit
import matplotlib.pyplot as plt
import scipy.linalg as scla
import matplotlib.gridspec as gridspec
import scipy.interpolate
import os

# Parameters for object a
flux_a = 0            # We will loop through this value later
hlr_a = 1             # arcsec
e1_a = 0.0
e2_a = 0.0
x0_a = -1
y0_a = 0
n_a = 0.5

# Parameters for object b
flux_b = 0            # We will loop through this value later
hlr_b = hlr_a         # arcsec
e1_b = 0.0
e2_b = 0.0
x0_b = 1
y0_b = 0
n_b = 0.5

true_val = np.array([flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b])

# Track the separation
sep = x0_b - x0_a

# Galsim function definitions
sersic_func = galsim.Sersic

# Set the RNG
seed_1 = galsim.BaseDeviate(0)
seed_2 = galsim.BaseDeviate(0)
seed_3 = galsim.BaseDeviate(0)

# Image properties
pixel_scale = 1/5     # arcsec / pixel
x_len = y_len = 100            # pixel

# Use LSST defined sky noise for r-band
add_noise_flag = True
texp = 6900 # seconds
sbar = 26.8 # sky photons per second per pixel
sky_level = texp*sbar
sky_noise = np.sqrt(sky_level)

# psf properties
psf_flag = False
beta = 3
fwhm_psf = 0.6

# Calc function for SNR --> flux
SNR_to_flux, snr_points, flux_pts = noiseLibrary.calc_SNR_to_flux(hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,
                                                        hlr_b,e1_b,e2_b,x0_b,y0_b,n_b,
                                                        False,beta,fwhm_psf,
                                                        x_len,y_len,pixel_scale,sersic_func,sersic_func,seed_1,seed_2,seed_3,
                                                        False,sky_level,sbar,texp,
                                                        10000,1000,3)
                                            
#plt.scatter(snr_points,flux_pts,c='b',alpha=0.5)
#plt.title('Flux vs SNR'); plt.xlabel('SNR'); plt.ylabel('Flux')
snr_points = np.array(snr_points); flux_pts = np.array(flux_pts) 
cond = np.logical_and(snr_points > 0, snr_points < 150)
flux_pts = flux_pts[cond]
snr_points = snr_points[cond]
#plt.xlim([0,np.max(snr_points)]); plt.ylim([0,np.max(flux_pts)])
SNR_to_flux = scipy.interpolate.interp1d(snr_points,flux_pts,kind='cubic')


# SNR range to loop through                                            
SNR_range = [100,40,30,20,15,10,5]
SNR = SNR_range[0]
# number of trials
num_trials = 500

# Data to keep track of
resid_matrix = []
error_matrix = []
pull_matrix = []

tot_flux = SNR_to_flux(SNR)
flux_a = tot_flux/2
flux_b = tot_flux/2
true_val[0] = flux_a; true_val[6] = flux_b

for i in xrange(0,num_trials):
    # Obtain instantiation
    im_no_noise, im_noise, best_fit, result = noiseLibrary.run_2_galaxy_full_params_simple(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,
                                                                                           flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,n_b,
                                                                                           psf_flag,beta,fwhm_psf,
                                                                                           x_len,y_len,pixel_scale,sersic_func,sersic_func,seed_1,seed_2,seed_3,
                                                                                           add_noise_flag,sky_level)
                                                                                           
    # Obtain the covariance and correlation matrix.
    if result.covar is None:
        print "Covariance Matrix is Null"
        im_no_noise, im_noise, best_fit, result = noiseLibrary.run_2_galaxy_full_params_simple(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,
                                                                                               flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,n_b,
                                                                                               psf_flag,beta,fwhm_psf,
                                                                                               x_len,y_len,pixel_scale,sersic_func,sersic_func,seed_1,seed_2,seed_3,
                                                                                               add_noise_flag,sky_level)
                                                                                   
    error_diag = (np.sqrt(np.diag(result.covar)))
    error_mat = np.outer(error_diag,error_diag)
    correlation_mat = result.covar/error_mat
    error_diag = error_diag.tolist()
    
    # Obtain the residuals
    resid = noiseLibrary.calc_resid(result.params,true_val)
    
    # Put residuals and error bars in a matrix and obtain the pull
    pull = np.array(resid)/error_diag
    # Put pull in this form in order to plot with imshow
    pull_mat = np.array([[pull[0]],[pull[1]],[pull[2]],[pull[3]],[pull[4]],[pull[5]],
                         [pull[6]],[pull[7]],[pull[8]],[pull[9]],[pull[10]],[pull[11]]])                                                                                       
    pull = pull.tolist()
    # Store the data in vectors
    resid_matrix.append(resid)
    error_matrix.append(error_diag)
    pull_matrix.append(pull)

resid_matrix = np.array(resid_matrix)
error_matrix = np.array(error_matrix)
pull_matrix = np.array(pull_matrix)

file = np.savetxt('SNR_' + str(SNR) + '_sep_' + str(sep) + '_trials_' + str(num_trials) + '.txt',resid_matrix)

