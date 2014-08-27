# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 13:54:38 2014

@author: luis
"""

## This script runs through a simple double fit of two sersic profiles

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
flux_a = 100000*0.9          # total counts on the image
hlr_a = 1             # arcsec
e1_a = 0.0
e2_a = 0.0
x0_a = -1
y0_a = 0
n_a = 0.5

# Parameters for object b
flux_b = flux_a          # total counts on the image
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
seed_1 = galsim.BaseDeviate(1)
seed_2 = galsim.BaseDeviate(2)
seed_3 = galsim.BaseDeviate(3)

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

# Obtain instantiation
im_no_noise, im_noise, best_fit, result = noiseLibrary.run_2_galaxy_full_params_simple(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,
                                                                                       flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,n_b,
                                                                                       psf_flag,beta,fwhm_psf,
                                                                                       x_len,y_len,pixel_scale,sersic_func,sersic_func,seed_1,seed_2,seed_3,
                                                                                       add_noise_flag,sky_level)

# Calc SNR for the noise-free image.
nu_s, mask_per_s, nu, mask = noiseLibrary.calc_SNR(im_no_noise, texp, sbar, 0.5)

# Calc function for SNR --> flux
SNR_to_flux, snr_points, flux_pts = noiseLibrary.calc_SNR_to_flux(hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,
                                                        hlr_b,e1_b,e2_b,x0_b,y0_b,n_b,
                                                        psf_flag,beta,fwhm_psf,
                                                        x_len,y_len,pixel_scale,sersic_func,sersic_func,seed_1,seed_2,seed_3,
                                                        False,sky_level,sbar,texp,
                                                        10000,1000,10)
                                            
plt.figure()
plt.scatter(snr_points,flux_pts,c='g')
snr_points = np.array(snr_points); flux_pts = np.array(flux_pts) 
cond = np.logical_and(snr_points > 0, snr_points < 15)
flux_pts = flux_pts[cond]
snr_points = snr_points[cond]
plt.xlim([0,15]); plt.ylim([0,1e5])
SNR_to_flux = scipy.interpolate.interp1d(snr_points,flux_pts,kind='cubic')
plt.plot(snr_points,SNR_to_flux(snr_points),c='g')                                             
                                                                    
# Report errors
lmfit.report_errors(result.params)

# Obtain the covariance and correlation matrix.
error_diag = np.sqrt(np.diag(result.covar))
error_mat = np.outer(error_diag,error_diag)
correlation_mat = result.covar/error_mat

# Obtain the residuals
resid_a = noiseLibrary.calc_resid(result.params,true_val)

# Put residuals and error bars in a matrix and obtain the pull
pull = np.array(resid_a)/error_diag
# Put pull in this form in order to plot with imshow
pull = np.array([[pull[0]],[pull[1]],[pull[2]],[pull[3]],[pull[4]],[pull[5]],
                 [pull[6]],[pull[7]],[pull[8]],[pull[9]],[pull[10]],[pull[11]]])
                 
# Provide the fontsize
fonts = 10
# Initialize the geometry of the grid
gs = gridspec.GridSpec(3,12)
# Set the figure size
fig = plt.figure(figsize=(15,10))
plt.suptitle(r'$Image\/Scale:$ $\frac{0.2\rq\rq}{Pixel}$',fontsize=fonts+10)
# Add the first subplot
ax1 = fig.add_subplot(gs[0,:2])
plt.title('Image',fontsize=fonts+8)
plt.ylabel('$Pixels$',fontsize=fonts+8)
plt.tick_params(axis='both',which='major',labelsize=9)
a = ax1.imshow(im_noise.array,interpolation='none',origin='lower',vmax=im_noise.array.max(),vmin=im_noise.array.min())
plt.colorbar(a,shrink=1)
# Add the second subplot
ax2 = fig.add_subplot(gs[1,:2])
plt.ylabel('$Pixels$',fontsize=fonts+8)
plt.tick_params(axis='both',which='major',labelsize=9)
plt.title('Best Fit',fontsize=fonts+8)
b = ax2.imshow(best_fit.array,origin='lower',vmax=best_fit.array.max(),vmin=best_fit.array.min())
plt.colorbar(b,shrink=1)
# Add the third subplot
ax3 = fig.add_subplot(gs[2,:2])
plt.tick_params(axis='both',which='major',labelsize=9)
plt.title('Residual',fontsize=fonts+8)
plt.xlabel('$Pixels$',fontsize=fonts+8)
plt.ylabel('$Pixels$',fontsize=fonts+8)
c = ax3.imshow((best_fit - im_noise).array,interpolation='none',origin='lower')
plt.colorbar(c,shrink=1)
# Add the column of pulls
ax4 = fig.add_subplot(gs[:,3:4])
plt.title(r'Residual/Error',fontsize=fonts+6)
plt.yticks([0,1,2,3,4,5,6,7,8,9,10,11],['$Flux_a$','$HLR_a$','$e1_a$','$e2_a$','$x0_a$','$y0_a$',
           '$Flux_b$','$HLR_b$','$e1_b$','$e2_b$','$x0_b$','$y0_b$'],fontsize=15)
plt.xticks([0],['']) 
e = ax4.imshow(np.around(pull,decimals=2),interpolation='none',origin='lower',vmin=-3,vmax=3)
plt.colorbar(e,shrink=0.9)          
for (i,j), val in np.ndenumerate(pull):
    ax4.annotate('%0.2f'%(val), (j,i), ha='center', va='center',size=12)          
# Add the correlation matrix
ax5 = fig.add_subplot(gs[:,5:])
plt.title('Correlation Coefficient Matrix',fontsize=fonts+8)
plt.xticks([0,1,2,3,4,5,6,7,8,9,10,11],['$Flux_a$','$HLR_a$','$e1_a$','$e2_a$','$x0_a$','$y0_a$',
           '$Flux_b$','$HLR_b$','$e1_b$','$e2_b$','$x0_b$','$y0_b$'],fontsize=15)
plt.yticks([0,1,2,3,4,5,6,7,8,9,10,11],['$Flux_a$','$HLR_a$','$e1_a$','$e2_a$','$x0_a$','$y0_a$',
           '$Flux_b$','$HLR_b$','$e1_b$','$e2_b$','$x0_b$','$y0_b$'],fontsize=15)
d = ax5.imshow(np.around(correlation_mat,decimals=2),interpolation='none',origin='lower',vmin=-1,vmax=1)
plt.colorbar(d,shrink=0.8)
for (i,j), val in np.ndenumerate(correlation_mat):
    ax5.annotate('%0.2f'%(val), (j,i), ha='center', va='center',size=12)
text = ax4.text(0,-1.5,'$Separation\/=\/%.1f$\"; $Mixture\/SNR=\/%i$; $PSF_{fwhm}\/=\/0.6$\"'%(sep,nu),fontsize=15)

path = '/home/luis/Documents/SRC2014/galsim_work/Noise_Studies'
filename = 'fig.png'
filename = os.path.join(path, filename)
fig.savefig(filename)